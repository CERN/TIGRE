# Wrapper for Tigre

This folder contains:

* A command line interface to reconstruct tomograms and perform denoising with Tigre.
* Unitary tests end-to-end to warranty the tigre reconstruction

# Using Tigre through the command line

All Tigre algorithms can be executed via command line. However the command line wrapper is only valide for laminography. It considers a parallel beam, and the source and the detector are parallel, being the sample placed between them.a The alignment is given by two files: angles.tlt and inPlaneRotation.xf. The files .tlt contains the tilt angles of the sample. The file .xf contains the in plane rotations of the tilt images as well as the shifts. This two files follow the IMOD format of cryoET. The tilt images should be introduced in format .mrc

```
python3 tigre_reconstruction.py --tiltseries images.mrc --angles angles.tlt --xf inplaneRotations.xf --thickness 100 -o tomogram.mrc --method wbp --gpu 0
python3 tigre_reconstruction.py --tiltseries images.mrc --angles angles.tlt --xf inplaneRotations.xf --thickness 100 -o tomogram.mrc --method sirt --iter 20 --gpu 0
```


# Unitary test

All Tigre reconstruction algorithm have a tests end-to-end. The tests are located in the file `tests/test_cli.py`.
To execute the test just launch the next command line

```
python3 -m pytest test_cli.py
```
