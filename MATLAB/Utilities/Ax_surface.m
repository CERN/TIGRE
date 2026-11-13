function projections=Ax_surface(image,geo,surface,device)
%AX_SURFACE Native Siddon projection onto explicit curved/modular pixels.
% image: single/double [nx ny nz]; result: single [rows columns views].
% Coordinates are world xyz in mm. Masked gaps return zero. device is a
% zero-based CUDA device index (default 0); no gpuArray or Python is used.
if nargin<4, device=0; end
projections=surfaceProject(false,image,geo,surface,device);
end
