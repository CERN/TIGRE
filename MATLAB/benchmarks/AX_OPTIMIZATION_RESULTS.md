# Metal `Ax` optimization results

Measurements were collected on 10 September 2026 with MATLAB R2026a on an
Apple M2 Pro (16 GPU cores), macOS 15.5. The benchmark used two warmups and
five synchronized samples per case. Timings are medians.

## Retained changes

1. Siddon traversal now maintains a 32-bit linear voxel index and updates it
   when crossing voxel boundaries. On full d13 Siddon this changed the median
   from 0.4085 s to 0.3830 s (1.07x), with bit-identical output.
2. Page-aligned MATLAB input and output memory is wrapped directly in shared
   Metal buffers. A checked copy fallback remains for unaligned addresses or
   lengths. Combined with incremental indexing, d13 Siddon reached 0.3246 s
   (1.26x overall), with bit-identical output.
3. The same checked no-copy path is used by `Atb`. A standalone d13 matched
   backprojection changed from 0.2892 s to 0.2540 s (1.14x).

Opt-in instrumentation measured warmed d13 Siddon at approximately 300-311 ms
of Metal GPU execution after optimization. The no-copy path reduced image
setup from about 3 ms to 0.1 ms and eliminated a 52-55 ms output copy.
Geometry preparation remained below 0.1 ms.

## End-to-end benchmark

| Case | Projector | Before (s) | After (s) | Speedup |
|---|---:|---:|---:|---:|
| d13 helical | Siddon | 0.408489 | 0.324571 | 1.259x |
| d13 helical | interpolated | 2.958149 | 2.881596 | 1.027x |
| circular cone | Siddon | 0.011587 | 0.010431 | 1.111x |
| circular cone | interpolated | 0.090000 | 0.089644 | 1.004x |
| parallel | Siddon | 0.008728 | 0.006772 | 1.289x |
| parallel | interpolated | 0.061632 | 0.060019 | 1.027x |
| static-detector DBT | Siddon | 0.002901 | 0.002608 | 1.112x |
| static-detector DBT | interpolated | 0.013777 | 0.012246 | 1.125x |
| small overhead | Siddon | 0.002283 | 0.001352 | 1.689x |
| small overhead | interpolated | 0.002724 | 0.001359 | 2.005x |
| d13 CGLS, 20 requested iterations | — | 28.591648 | 21.443071 | 1.333x |

Every optimized benchmark output, including the d13 CGLS reconstruction, was
bit-for-bit identical to its frozen pre-optimization Metal output. The final
d13 Siddon projection had relative L2 error `3.228e-4` against CUDA. The
unchanged full d13 demo had CUDA relative L2 errors of `2.411e-4` for OS-SART
and `5.174e-4` for CGLS.

Baseline artifact:
`metal_ax_baseline_20260910T112936Z-8b8de33.mat`
(`47cc342cf68df63416d844a35de89741caea77bf25afea35c9414766500c9ff8`).

Optimized artifact:
`metal_ax_optimized_20260910T112936Z-8b8de33.mat`
(`f83c70f1297ec6ea01fe4b5b412197548af4023328b30be6b263669ec6630513`).

## Rejected experiments

- Detector-coordinate 3-D dispatch with `8x8`, `16x4`, and `16x8`
  threadgroups was bit-identical but did not outperform the existing
  one-dimensional dispatch.
- One-dimensional threadgroup widths from 32 through 1024 had overlapping
  timings and no repeatable improvement over 256.
- An exact-read `R32Float` 3-D texture path was bit-identical but regressed
  d13 Siddon to 0.5905 s, including texture creation and upload.
- Scratch-buffer reuse changed the tiny-case median by less than 5% and was
  not retained.

The Xcode Metal System Trace captured the command-buffer stream, but the stock
template had shader timeline and GPU counter collection disabled. Occupancy and
cache counters were therefore unavailable; command-buffer GPU timestamps and
end-to-end measurements drove the retained changes.
