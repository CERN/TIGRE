function [VM_proj] = MakeVMproj(Al_proj, PMMA_proj, energy)
%MAKEVMPROJ Synethesizes VM projections from basis material thickness
%projections.
%   
%   Coded by: Andrew Keeler
%   Contact: akeeler@luc.edu
%   Date: February 28, 2024

if (energy < 20 || energy > 150)
    warning('Requested energy is outside of expected range.  Results may be unreliable.')

%  linear attenuation values for basis materials derived from NIST XCOM
%  database
energies = [20, 30, 40, 50, 60, 80, 10, 110, 120, 130, 140, 150];
al_array = [0.9291 0.3046 0.1535 0.0994 0.0750 0.0545 0.0460 0.0431 0.0414 0.0398 0.0384 0.0372];
pmma_array = [0.0677 0.0359 0.0278 0.0246 0.0228 0.0207 0.0194 0.0188 0.0184 0.0179 0.0175 0.0172];

%  Compute attenuation coefficients for basis materials at selected energy
al_atten = spline(energies, al_array, energy);
pmma_atten = spline(energies, pmma_array, energy);

VM_proj = al_atten*Al_proj + pmma_atten*PMMA_proj;
end

