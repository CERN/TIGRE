# TIGRE 3.3.1 vs ASTRA: speed of the projection operators

`tigre331_vs_astra.pdf` is the technical note; `tigre331_vs_astra.tex` its source, with its tables in `tables/`.

## Run the comparison on your machine

You need TIGRE (3.3.1 or later), the ASTRA Toolbox (2.5), PyTorch with CUDA, and NumPy, in one Python
environment. Keep `tigre_vs_astra_suite.py` and `tigre_vs_astra.py` in the same folder, and run from outside
TIGRE's source tree:

```
python tigre_vs_astra_suite.py --check                      # find ASTRA's conventions; quick
python tigre_vs_astra_suite.py --ops fp bp --gpus 1 --label my-gpu --out suite.json --csv suite.csv
```

- `--gpus 1 2 4` also measures several GPUs (as many as you have).
- `--groups paper --ops fp bp --gpus 1 2 4` runs the large problems of the note: N^3 images and N angles, up to
  1536^3, which needs about 70 GB of RAM.
- `--list` prints the cases without running them, and `--resume` carries on after a stop.

The suite takes about an hour on one GPU, and several hours with `--gpus 1 2 4`.

**Results from other GPUs and multi-GPU machines are very welcome.** The note's multi-GPU results come from a
single server, and the thresholds TIGRE uses (Siddon steps per pass, #916; page-locking from three GPUs, #917)
come from two GPU models. Please open an issue on GitHub with your `suite.json` and a line about the machine:
GPUs, how they are connected (`nvidia-smi topo -m`), CPU, RAM and OS.

## Remake the tables

`results/` holds the raw results of the note: `win_*` from the RTX 4070 (Windows) and `linux_*` from the
Quadro RTX 6000 server, for 3.3.1 and for 3.3.0. (The 3.3.1 runs were made on master with #911–#917 merged, before the
version number was raised, so their `meta.tigre` still reads 3.3.0.) `./regen.sh` rebuilds `tables/` from them.
`python3 losses.py results/linux_331.json 1` lists every comparison that ASTRA wins on one GPU.
