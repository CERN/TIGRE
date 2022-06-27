import matplotlib.patches
import numpy as np
import tigre
from mpl_toolkits.mplot3d import art3d

# https://matplotlib.org/mpl_toolkits/mplot3d/tutorial.html


class Arrow3D(matplotlib.patches.FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        matplotlib.patches.FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        from mpl_toolkits.mplot3d import proj3d

        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        matplotlib.patches.FancyArrowPatch.draw(self, renderer)

    def do_3d_projection(self, render=None):
        from mpl_toolkits.mplot3d import proj3d
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return np.min(zs)


ROT_DEFAULT = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
TRANS_DEFAULT = np.array([0, 0, 0])


def pathpatch_2d_to_3d_affine(pathpatch, mat_rot=ROT_DEFAULT, vec_trans=TRANS_DEFAULT):
    """
    Transforms a 2D Patch to a 3D patch using the affine tranform
    of the given rotation matrix and translation vector.
    The pathpatch is assumed to be on the plane Z = 0.
    """
    path = pathpatch.get_path()  # Get the path and the associated transform
    trans = pathpatch.get_patch_transform()

    path = trans.transform_path(path)  # Apply the transform

    pathpatch.__class__ = art3d.PathPatch3D  # Change the class
    pathpatch._code3d = path.codes  # Copy the codes
    pathpatch._facecolor3d = pathpatch.get_facecolor  # Get the face color

    verts = path.vertices  # Get the vertices in 2D

    M = np.array(
        [
            [mat_rot[0, 0], mat_rot[0, 1], mat_rot[0, 2], vec_trans[0]],
            [mat_rot[1, 0], mat_rot[1, 1], mat_rot[1, 2], vec_trans[1]],
            [mat_rot[2, 0], mat_rot[2, 1], mat_rot[2, 2], vec_trans[2]],
        ]
    )

    pathpatch._segment3d = np.array([np.dot(M, (x, y, 0, 1)) for x, y in verts])


def plot_geometry(geo, angle=0):
    """Plots the given geometry."""
    import mpl_toolkits.mplot3d.art3d as art3d
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(6, 6))
    fig.suptitle("Cone Beam Compute Tomography geometry")
    ax = fig.add_subplot(111, projection="3d")
    ax.set_title("Current CBCT geometry, in scale")

    limXY = max(geo.DSO, geo.DSD - geo.DSO)
    limZ = geo.sVoxel[0]

    ax.set_box_aspect((limXY, limXY, limZ))

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_xlim3d(-limXY * 1.2, limXY * 1.2)
    ax.set_ylim3d(-limXY * 1.2, limXY * 1.2)
    ax.set_zlim3d(-limZ * 1.2, limZ * 1.2)

    # Trajectory of Source
    # https://matplotlib.org/devdocs/api/_as_gen/matplotlib.patches.Circle.htm
    circ = matplotlib.patches.Circle((0, 0), geo.DSO, color="black", fill=False, ls="-.", lw=0.5)
    ax.add_patch(circ)
    art3d.pathpatch_2d_to_3d(circ, z=0, zdir="z")

    # Trajectory of Detector
    circ = matplotlib.patches.Circle(
        (0, 0), geo.DSD - geo.DSO, color="black", fill=False, ls="-.", lw=0.5
    )
    ax.add_patch(circ)
    art3d.pathpatch_2d_to_3d(circ, z=0, zdir="z")

    # Source
    # ax.scatter([0], [0], [0], color="g", s=100)
    sourcePos3D = [geo.DSO * np.cos(angle), geo.DSO * np.sin(angle), 0]  # xyz
    ax.scatter([sourcePos3D[0]], [sourcePos3D[1]], [sourcePos3D[2]], color="steelblue", s=100)

    # Axes XYZ
    length_axis = geo.sVoxel[0]
    x_axis = Arrow3D(
        [0, length_axis],
        [0, 0],
        [0, 0],
        mutation_scale=10,
        shrinkA=0,
        lw=1,
        arrowstyle="-|>",
        color="r",
    )
    y_axis = Arrow3D(
        [0, 0],
        [0, length_axis],
        [0, 0],
        mutation_scale=10,
        shrinkA=0,
        lw=1,
        arrowstyle="-|>",
        color="b",
    )
    z_axis = Arrow3D(
        [0, 0],
        [0, 0],
        [0, length_axis],
        mutation_scale=10,
        shrinkA=0,
        lw=1,
        arrowstyle="-|>",
        color="g",
    )
    ax.add_artist(x_axis)
    ax.add_artist(y_axis)
    ax.add_artist(z_axis)

    # Detector
    print("sDetector = {}".format(geo.sDetector))
    alpha_detector = 0.7
    detectorPos3D = [
        (geo.DSD - geo.DSO) * np.cos(angle + np.pi) + geo.offDetector[1] * np.sin(angle + np.pi),
        (geo.DSD - geo.DSO) * np.sin(angle + np.pi) - geo.offDetector[1] * np.cos(angle + np.pi),
        geo.offDetector[0],
    ]  # xyz
    detector = matplotlib.patches.Rectangle(
        (-geo.sDetector[1] / 2, -geo.sDetector[0] / 2),
        geo.sDetector[1],
        geo.sDetector[0],
        angle=0,  # angle,
        color="maroon",
        fill=True,
        alpha=alpha_detector,
    )
    print("detector={}".format(detector))
    ax.add_patch(detector)
    mat_rot = np.array(
        [[-np.sin(angle), 0, np.cos(angle)], [np.cos(angle), 0, np.sin(angle)], [0, 1, 0]]
    )
    pathpatch_2d_to_3d_affine(detector, mat_rot, detectorPos3D)

    # Image FOV
    alpha_img = 0.1
    offOrigin = np.array([geo.offOrigin[2], geo.offOrigin[1], geo.offOrigin[0]])
    mat_rot_xy = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float)
    for idx in range(2):
        img_face_xy = matplotlib.patches.Rectangle(
            (-geo.sVoxel[2] / 2, -geo.sVoxel[1] / 2),  # xy order
            geo.sVoxel[2],
            geo.sVoxel[1],
            angle=0,  # angle,
            color="black",
            fill=True,
            alpha=alpha_img,
        )
        ax.add_patch(img_face_xy)
        pathpatch_2d_to_3d_affine(
            img_face_xy,
            mat_rot_xy,
            np.array([0, 0, -geo.sVoxel[0] / 2 if idx == 0 else geo.sVoxel[0] / 2]) + offOrigin,
        )

    mat_rot_yz = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=np.float)
    for idx in range(2):
        img_face_yz = matplotlib.patches.Rectangle(
            (-geo.sVoxel[1] / 2, -geo.sVoxel[0] / 2),  # xy order
            geo.sVoxel[1],
            geo.sVoxel[0],
            angle=0,  # angle,
            color="black",
            fill=True,
            alpha=alpha_img,
        )
        ax.add_patch(img_face_yz)
        pathpatch_2d_to_3d_affine(
            img_face_yz,
            mat_rot_yz,
            np.array([-geo.sVoxel[2] / 2 if idx == 0 else geo.sVoxel[2] / 2, 0, 0]) + offOrigin,
        )

    mat_rot_zx = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=np.float)
    for idx in range(2):
        img_face_zx = matplotlib.patches.Rectangle(
            (-geo.sVoxel[0] / 2, -geo.sVoxel[2] / 2),  # xy order
            geo.sVoxel[0],
            geo.sVoxel[2],
            angle=0,  # angle,
            color="black",
            fill=True,
            alpha=alpha_img,
        )
        ax.add_patch(img_face_zx)
        pathpatch_2d_to_3d_affine(
            img_face_zx,
            mat_rot_zx,
            np.array([0, -geo.sVoxel[1] / 2 if idx == 0 else geo.sVoxel[1] / 2, 0]) + offOrigin,
        )

    plt.show()


if __name__ == "__main__":
    n_voxel_z = 128
    n_voxel_y = 256
    n_voxel_x = 512
    n_detector_u = 400
    n_detector_v = 300
    off_detector_u = 0  # 500
    off_detector_v = 0  # 500
    off_origin_x = 0  # 300
    off_origin_y = 0  # 100
    off_origin_z = 0  # 100

    geo = tigre.geometry(mode="cone", default=True)
    geo.nVoxel = np.array([n_voxel_z, n_voxel_y, n_voxel_x])
    geo.sVoxel = geo.nVoxel
    geo.dVoxel = geo.sVoxel / geo.nVoxel
    geo.nDetector = np.array([n_detector_v, n_detector_u])
    geo.sDetector = geo.nDetector * geo.dDetector
    geo.offDetector = np.array([off_detector_v, off_detector_u])
    geo.offOrigin = np.array([off_origin_z, off_origin_y, off_origin_x])
    print(geo)
    angle = -np.pi / 6
    plot_geometry(geo, angle)
