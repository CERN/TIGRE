function image=Atb_surface(projections,geo,surface,device)
%ATB_SURFACE Discrete adjoint of Ax_surface, including its valid_mask.
% projections: [rows columns views]; result: single [nx ny nz].
% This is a matched backprojector, not an FDK reconstruction.
% Float32 atomic accumulation is subject to summation-order rounding. Masked
% NaNs are discarded; nonfinite valid samples are rejected. Input is unchanged.
if nargin<4, device=0; end
image=surfaceProject(true,projections,geo,surface,device);
end
