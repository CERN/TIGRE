from tigre.algorithms import OSSART
from tigre.utilities import example_geometry, example_projections

geo = example_geometry()
proj = example_projections()

# Run without Parker weighting
vol1 = OSSART(proj, geo, angles=None, niter=1)

# Run with Parker weighting
vol2 = OSSART(proj, geo, angles=None, niter=1, parker_weighting=True)

print("OSSART ran successfully with and without Parker weighting!")
