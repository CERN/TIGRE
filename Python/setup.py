import copy
import glob
import os
from os.path import join as pjoin
import re
import subprocess
import sys
import shutil
import re

from Cython.Distutils import build_ext
import numpy
from setuptools import setup, find_packages, Extension


IS_WINDOWS = sys.platform == "win32"


# Code from https://github.com/pytorch/pytorch/blob/master/torch/utils/cpp_extension.py
COMPUTE_CAPABILITY_ARGS = [  # '-gencode=arch=compute_20,code=sm_20', #deprecated
    "-gencode=arch=compute_30,code=sm_30",  # Deprecated at CUDA 9.2
    "-gencode=arch=compute_37,code=sm_37",
    "-gencode=arch=compute_50,code=sm_50",
    "-gencode=arch=compute_52,code=sm_52",
    "-gencode=arch=compute_60,code=sm_60",
    "-gencode=arch=compute_61,code=sm_61",
    "-gencode=arch=compute_70,code=sm_70",
    "-gencode=arch=compute_75,code=sm_75",  # From CUDA 10
    "-gencode=arch=compute_86,code=sm_86",  # From CUDA 11
    "-gencode=arch=compute_70,code=compute_70", # allows foward compiling
    "--ptxas-options=-v",
    "-c",
    "--default-stream=per-thread",
]


def get_cuda_version(cuda_home):
    """Locate the CUDA version"""
    version_file = os.path.join(cuda_home, "version.txt")
    try:
        if os.path.isfile(version_file):
            with open(version_file) as f:
                version_str = f.readline().replace("\n", "").replace("\r", "")
                return version_str.split(" ")[2][:4]
        else:
            version_str = subprocess.check_output(
                [os.path.join(cuda_home, "bin", "nvcc"), "--version"]
            )
            version_str = str(version_str).replace("\n", "").replace("\r", "")
            idx = version_str.find("release")
            return version_str[idx + len("release ") : idx + len("release ") + 4]
    except:
        raise RuntimeError("Cannot read cuda version file")


def locate_cuda():
    """Locate the CUDA environment on the system

    Returns a dict with keys 'home', 'include' and 'lib64'
    and values giving the absolute path to each directory.

    Starts by looking for the CUDA_HOME or CUDA_PATH env variable. If not found, everything
    is based on finding 'nvcc' in the PATH.
    """
    # Guess #1
    cuda_home = os.environ.get("CUDA_HOME") or os.environ.get("CUDA_PATH")
    if cuda_home is None:
        # Guess #2
        try:
            which = "where" if IS_WINDOWS else "which"
            nvcc = subprocess.check_output([which, "nvcc"]).decode().rstrip("\r\n")
            cuda_home = os.path.dirname(os.path.dirname(nvcc))
        except subprocess.CalledProcessError:
            # Guess #3
            if IS_WINDOWS:
                cuda_homes = glob.glob("C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v*.*")
                if len(cuda_homes) == 0:
                    cuda_home = ""
                else:
                    cuda_home = cuda_homes[0]
            else:
                cuda_home = "/usr/local/cuda"
            if not os.path.exists(cuda_home):
                cuda_home = None
    version = get_cuda_version(cuda_home)
    cudaconfig = {
        "home": cuda_home,
        "include": pjoin(cuda_home, "include"),
        "lib64": pjoin(cuda_home, pjoin("lib", "x64") if IS_WINDOWS else "lib64"),
    }
    if not all([os.path.exists(v) for v in cudaconfig.values()]):
        raise EnvironmentError(
            "The CUDA  path could not be located in $PATH, $CUDA_HOME or $CUDA_PATH. "
            "Either add it to your path, or set $CUDA_HOME or $CUDA_PATH."
        )

    return cudaconfig, version


def _is_cuda_file(path):
    return os.path.splitext(path)[1] in [".cu", ".cuh"]


CUDA, CUDA_VERSION = locate_cuda()

try:
    cuda_version = float(CUDA_VERSION)
except ValueError:
    cuda_list = re.findall('\d+', CUDA_VERSION)
    cuda_version = float( str(cuda_list[0] + '.' + cuda_list[1]))

# Cleanup CUDA arguments depedning on the version
if cuda_version < 11.0:
    COMPUTE_CAPABILITY_ARGS.pop(8)

if cuda_version < 10.0:
    COMPUTE_CAPABILITY_ARGS.pop(7)

if cuda_version > 9.2:
    COMPUTE_CAPABILITY_ARGS.pop(0)
    
# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    NUMPY_INCLUDE = numpy.get_include()
except AttributeError:
    NUMPY_INCLUDE = numpy.get_numpy_include()


COMMON_MSVC_FLAGS = ["/MD", "/wd4819", "/EHsc"]


COMMON_NVCC_FLAGS = [
    "-D__CUDA_NO_HALF_OPERATORS__",
    "-D__CUDA_NO_HALF_CONVERSIONS__",
    "-D__CUDA_NO_HALF2_OPERATORS__",
    "--expt-relaxed-constexpr",
]


def _join_cuda_home(*paths):
    return os.path.join(CUDA["home"], *paths)


class BuildExtension(build_ext):
    """
    A custom :mod:`Cython.Distutils` build extension .

    This :class:`Cython.Distutils.build_ext` subclass takes care of passing the
    minimum required compiler flags (e.g. ``-std=c++11``) as well as mixed
    C++/CUDA compilation (and support for CUDA files in general).

    When using :class:`BuildExtension`, it is allowed to supply a dictionary
    for ``extra_compile_args`` (rather than the usual list) that maps from
    languages (``cxx`` or ``nvcc``) to a list of additional compiler flags to
    supply to the compiler. This makes it possible to supply different flags to
    the C++ and CUDA compiler during mixed compilation.
    """

    @classmethod
    def with_options(cls, **options):
        """
        Returns an alternative constructor that extends any original keyword
        arguments to the original constructor with the given options.
        """

        def init_with_options(*args, **kwargs):
            kwargs = kwargs.copy()
            kwargs.update(options)
            return cls(*args, **kwargs)

        return init_with_options

    def __init__(self, *args, **kwargs):
        build_ext.__init__(self, *args, **kwargs)
        self.no_python_abi_suffix = kwargs.get("no_python_abi_suffix", False)

    def build_extensions(self):
        # Register .cu and .cuh as valid source extensions.
        self.compiler.src_extensions += [".cu", ".cuh"]
        # Save the original _compile method for later.
        if self.compiler.compiler_type == "msvc":
            self.compiler._cpp_extensions += [".cu", ".cuh"]
            original_compile = self.compiler.compile
            original_spawn = self.compiler.spawn
        else:
            original_compile = self.compiler._compile

        def unix_wrap_compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
            # Copy before we make any modifications.
            cflags = copy.deepcopy(extra_postargs)
            try:
                original_compiler = self.compiler.compiler_so
                if _is_cuda_file(src):
                    nvcc = _join_cuda_home("bin", "nvcc")
                    if not isinstance(nvcc, list):
                        nvcc = [nvcc]
                    self.compiler.set_executable("compiler_so", nvcc)
                    if isinstance(cflags, dict):
                        cflags = cflags["nvcc"]
                    cflags = (
                        COMMON_NVCC_FLAGS
                        + ["--compiler-options", "'-fPIC'"]
                        + cflags
                        + COMPUTE_CAPABILITY_ARGS
                    )
                elif isinstance(cflags, dict):
                    cflags = cflags["cxx"]
                # NVCC does not allow multiple -std to be passed, so we avoid
                # overriding the option if the user explicitly passed it.
                if not any(flag.startswith("-std=") for flag in cflags):
                    cflags.append("-std=c++11")

                original_compile(obj, src, ext, cc_args, cflags, pp_opts)
            finally:
                # Put the original compiler back in place.
                self.compiler.set_executable("compiler_so", original_compiler)

        def win_wrap_compile(
            sources,
            output_dir=None,
            macros=None,
            include_dirs=None,
            debug=0,
            extra_preargs=None,
            extra_postargs=None,
            depends=None,
        ):

            cflags = copy.deepcopy(extra_postargs)
            extra_postargs = None

            def spawn(cmd, cflags):
                # Using regex to match src, obj and include files
                src_regex = re.compile("/T(p|c)(.*)")
                src_list = [m.group(2) for m in (src_regex.match(elem) for elem in cmd) if m]

                obj_regex = re.compile("/Fo(.*)")
                obj_list = [m.group(1) for m in (obj_regex.match(elem) for elem in cmd) if m]

                include_regex = re.compile(r"((\-|\/)I.*)")
                include_list = [
                    m.group(1) for m in (include_regex.match(elem) for elem in cmd) if m
                ]

                if len(src_list) >= 1 and len(obj_list) >= 1:
                    src = src_list[0]
                    obj = obj_list[0]
                    if _is_cuda_file(src):
                        nvcc = _join_cuda_home("bin", "nvcc")
                        if isinstance(cflags, dict):
                            cflags = cflags["nvcc"]
                        elif not isinstance(cflags, list):
                            cflags = []

                        cflags = COMMON_NVCC_FLAGS + cflags + COMPUTE_CAPABILITY_ARGS
                        for flag in COMMON_MSVC_FLAGS:
                            cflags = ["-Xcompiler", flag] + cflags
                        for macro in macros:
                            if len(macro) == 2:
                                if macro[1] == None:
                                    cflags += ["--define-macro", macro[0]]
                                else:
                                    cflags += ["--define-macro", "{}={}".format(macro[0], macro[1])]
                            elif len(macro) == 1:
                                cflags += ["--undefine-macro", macro[0]]

                        cmd = [nvcc, "-c", src, "-o", obj] + include_list + cflags
                    elif isinstance(cflags, dict):
                        cflags = COMMON_MSVC_FLAGS  # + self.cflags['cxx']
                        cmd += cflags
                    elif isinstance(cflags, list):
                        cflags = COMMON_MSVC_FLAGS + cflags
                        cmd += cflags

                return original_spawn(cmd)

            try:
                self.compiler.spawn = lambda cmd: spawn(cmd, cflags)
                return original_compile(
                    sources,
                    output_dir,
                    macros,
                    include_dirs,
                    debug,
                    extra_preargs,
                    extra_postargs,
                    depends,
                )
            finally:
                self.compiler.spawn = original_spawn

        # Monkey-patch the _compile method.
        if self.compiler.compiler_type == "msvc":
            self.compiler.compile = win_wrap_compile
        else:
            self.compiler._compile = unix_wrap_compile

        build_ext.build_extensions(self)

    def get_ext_filename(self, ext_name):
        # Get the original shared library name. For Python 3, this name will be
        # suffixed with "<SOABI>.so", where <SOABI> will be something like
        # cpython-37m-x86_64-linux-gnu. On Python 2, there is no such ABI name.
        # The final extension, .so, would be .lib/.dll on Windows of course.
        ext_filename = build_ext.get_ext_filename(self, ext_name)
        # If `no_python_abi_suffix` is `True`, we omit the Python 3 ABI
        # component. This makes building shared libraries with setuptools that
        # aren't Python modules nicer.
        if self.no_python_abi_suffix and sys.version_info >= (3, 0):
            # The parts will be e.g. ["my_extension", "cpython-37m-x86_64-linux-gnu", "so"].
            ext_filename_parts = ext_filename.split(".")
            # Omit the second to last element.
            without_abi = ext_filename_parts[:-2] + ext_filename_parts[-1:]
            ext_filename = ".".join(without_abi)
        return ext_filename


def include_headers(filename_list, sdist=False):
    """add hpp and h files to list if sdist is called"""
    if not sdist:
        return filename_list

    c_extensions = [".cu", ".c", ".C", ".cc", ".cpp", ".cxx", ".c++"]
    header_list = []
    for filename in filename_list:
        header = list(os.path.splitext(filename))
        if header[1] in c_extensions:
            header[1] = ".hpp"
            header_list.append("".join(header))

    filename_list += ["../Common/CUDA/types_TIGRE.hpp", "../Common/CUDA/errors.hpp"]
    return filename_list + header_list


Ax_ext = Extension(
    "_Ax",
    sources=include_headers(
        [
            "../Common/CUDA/projection.cpp",
            "../Common/CUDA/TIGRE_common.cpp",
            "../Common/CUDA/Siddon_projection.cu",
            "../Common/CUDA/Siddon_projection_parallel.cu",
            "../Common/CUDA/ray_interpolated_projection.cu",
            "../Common/CUDA/ray_interpolated_projection_parallel.cu",
            "../Common/CUDA/GpuIds.cpp",
            "tigre/utilities/cuda_interface/_types.pxd",
            "tigre/utilities/cuda_interface/_gpuUtils.pxd",
            "tigre/utilities/cuda_interface/_Ax.pyx",
        ],
        sdist=sys.argv[1] == "sdist",
    ),
    define_macros=[("IS_FOR_PYTIGRE", None)],
    library_dirs=[CUDA["lib64"]],
    libraries=["cudart"],
    language="c++",
    runtime_library_dirs=[CUDA["lib64"]] if not IS_WINDOWS else None,
    include_dirs=[NUMPY_INCLUDE, CUDA["include"], "../Common/CUDA/"],
)


Atb_ext = Extension(
    "_Atb",
    sources=include_headers(
        [
            "../Common/CUDA/TIGRE_common.cpp",
            "../Common/CUDA/voxel_backprojection.cu",
            "../Common/CUDA/voxel_backprojection2.cu",
            "../Common/CUDA/voxel_backprojection_parallel.cu",
            "../Common/CUDA/GpuIds.cpp",
            "../Common/CUDA/gpuUtils.cu",
            "tigre/utilities/cuda_interface/_types.pxd",
            "tigre/utilities/cuda_interface/_Atb.pyx",
        ],
        sdist=sys.argv[1] == "sdist",
    ),
    define_macros=[("IS_FOR_PYTIGRE", None)],
    library_dirs=[CUDA["lib64"]],
    libraries=["cudart"],
    language="c++",
    runtime_library_dirs=[CUDA["lib64"]] if not IS_WINDOWS else None,
    include_dirs=[NUMPY_INCLUDE, CUDA["include"], "../Common/CUDA/"],
)


tvdenoising_ext = Extension(
    "_tvdenoising",
    sources=include_headers(
        [
            "../Common/CUDA/TIGRE_common.cpp",
            "../Common/CUDA/tvdenoising.cu",
            "../Common/CUDA/GpuIds.cpp",
            "../Common/CUDA/gpuUtils.cu",
            "tigre/utilities/cuda_interface/_types.pxd",
            "tigre/utilities/cuda_interface/_tvdenoising.pyx",
        ],
        sdist=sys.argv[1] == "sdist",
    ),
    define_macros=[("IS_FOR_PYTIGRE", None)],
    library_dirs=[CUDA["lib64"]],
    libraries=["cudart"],
    language="c++",
    runtime_library_dirs=[CUDA["lib64"]] if not IS_WINDOWS else None,
    include_dirs=[NUMPY_INCLUDE, CUDA["include"], "../Common/CUDA/"],
)


minTV_ext = Extension(
    "_minTV",
    sources=include_headers(
        [
            "../Common/CUDA/TIGRE_common.cpp",
            "../Common/CUDA/POCS_TV.cu",
            "../Common/CUDA/GpuIds.cpp",
            "../Common/CUDA/gpuUtils.cu",
            "tigre/utilities/cuda_interface/_types.pxd",
            "tigre/utilities/cuda_interface/_minTV.pyx",
        ],
        sdist=sys.argv[1] == "sdist",
    ),
    define_macros=[("IS_FOR_PYTIGRE", None)],
    library_dirs=[CUDA["lib64"]],
    libraries=["cudart"],
    language="c++",
    runtime_library_dirs=[CUDA["lib64"]] if not IS_WINDOWS else None,
    include_dirs=[NUMPY_INCLUDE, CUDA["include"], "../Common/CUDA/"],
)


AwminTV_ext = Extension(
    "_AwminTV",
    sources=include_headers(
        [
            "../Common/CUDA/TIGRE_common.cpp",
            "../Common/CUDA/POCS_TV2.cu",
            "../Common/CUDA/GpuIds.cpp",
            "../Common/CUDA/gpuUtils.cu",
            "tigre/utilities/cuda_interface/_AwminTV.pyx",
        ],
        sdist=sys.argv[1] == "sdist",
    ),
    define_macros=[("IS_FOR_PYTIGRE", None)],
    library_dirs=[CUDA["lib64"]],
    libraries=["cudart"],
    language="c++",
    runtime_library_dirs=[CUDA["lib64"]] if not IS_WINDOWS else None,
    include_dirs=[NUMPY_INCLUDE, CUDA["include"], "../Common/CUDA/"],
)


gpuUtils_ext = Extension(
    "_gpuUtils",
    sources=include_headers(
        [
            "../Common/CUDA/gpuUtils.cu",
            "tigre/utilities/cuda_interface/_gpuUtils.pxd",
            "tigre/utilities/cuda_interface/_gpuUtils.pyx",
        ],
        sdist=sys.argv[1] == "sdist",
    ),
    library_dirs=[CUDA["lib64"]],
    libraries=["cudart"],
    language="c++",
    runtime_library_dirs=[CUDA["lib64"]] if not IS_WINDOWS else None,
    include_dirs=[NUMPY_INCLUDE, CUDA["include"], "../Common/CUDA/"],
)


RandomNumberGenerator_ext = Extension(
    "_RandomNumberGenerator",
    sources=include_headers(
        [
            "../Common/CUDA/TIGRE_common.cpp",
            "../Common/CUDA/RandomNumberGenerator.cu",
            "../Common/CUDA/GpuIds.cpp",
            "../Common/CUDA/gpuUtils.cu",
            "tigre/utilities/cuda_interface/_randomNumberGenerator.pyx",
        ],
        sdist=sys.argv[1] == "sdist",
    ),
    define_macros=[("IS_FOR_PYTIGRE", None)],
    library_dirs=[CUDA["lib64"]],
    libraries=["cudart"],
    language="c++",
    runtime_library_dirs=[CUDA["lib64"]] if not IS_WINDOWS else None,
    include_dirs=[NUMPY_INCLUDE, CUDA["include"], "../Common/CUDA/"],
)


setup(
    name="pytigre",
    version="2.2.0",
    author="Ander Biguri, Reuben Lindroos, Sam Loescher",
    packages=find_packages(),
    include_package_data=True,
    data_files=[("data", ["../Common/data/head.mat"])],
    ext_modules=[Ax_ext, Atb_ext, tvdenoising_ext, minTV_ext, AwminTV_ext, gpuUtils_ext, RandomNumberGenerator_ext],
    py_modules=["tigre.py"],
    cmdclass={"build_ext": BuildExtension},
    install_requires=["Cython", "matplotlib", "numpy", "scipy", "tqdm"],
    license_files=("LICENSE",),
    license="BSD 3-Clause",
    # since the package has c code, the egg cannot be zipped
    zip_safe=False,
)