# Python API

```python
import tigre                    # operators, geometry, plotting
import tigre.algorithms as algs # reconstruction algorithms
```

Arrays follow NumPy's index order: images are `(z, y, x)` and projections
`(angle, v, u)`, and the vector fields of the geometry are reversed compared
with MATLAB. All arrays passed to the GPU must be `float32`.

```{toctree}
:maxdepth: 1

operators
geometry
algorithms
utilities
io
plotting
```
