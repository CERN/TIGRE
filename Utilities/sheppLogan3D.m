function [p,ellipse]= sheppLogan3D( varargin )
%SHEPPLOGAN3D wrapper to phantom3DAniso, in order to maintain the original
%name, but change it in the toolbox
if isempty(varargin)
    [p,ellipse]=phantom3dAniso([128,128,128]);
else
    [p,ellipse]=phantom3dAniso(varargin);
end

end

