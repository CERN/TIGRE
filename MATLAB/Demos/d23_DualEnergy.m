%% Demo 23: Dual Energy CBCT in TIGRE. 
%
% This demo will illustrate the options for performing dual-energy analysis
% using the TIGRE toolkit.
% 
%   Coded by: Andrew Keeler
%   Contact: akeeler@luc.edu
%   Date: February 28, 2024

%% Data loading (See d21_ScannerDataLoader.m)
%   We'll use VarianDataLoader as an example.  Load two sets of scan data
%   at different kV settings into memory. Geo should be the same for both
%   scans, and so only needs to be loaded once.  Beam hardening correction
%   is not needed.

datafolder_low='~/your_data_path/varian/low_kV/2020-01-01_123456/';
datafolder_high='~/your_data_path/varian/high_kV/2020-01-01_123456/';
[proj_h, geo, angles_h] = VarianDataLoader(datafolder_high, 'bh', false);
[proj_l, ~, angles_l] = VarianDataLoader(datafolder, 'bh', false);

%% Material Decomposition
% DE-CBCT decomposes the two scans into equivalent thicknesses of two basis
% materials.  Generally, either Al and PMMA or Iodine and Water are used.
% The implementation here uses Al and PMMA.
%
% The decomposition procedure defaults to taking scans at 80 and 140 kV.
% Different energies can be specified using the 'kVl' and 'kVh' tags.

[Al_proj, PMMA_proj, angles] = DeDecompose(proj_l, angles_l, proj_h, angles_h);

%% Image construction
% Equivalent thickness projections can now be used to synthesize various
% DE-derived image types.  The most common of these are virtual
% monoenergetic (VM) images.

energy = 60; % keV
VM_proj = MakeVMproj(Al_proj, PMMA_proj, energy);

VMimg = FDK(VM_proj, geo, angles);
VMimg = OS_SART(VM_proj, geo, angles, 100);

%% Plot image
plotImg(VMimg)