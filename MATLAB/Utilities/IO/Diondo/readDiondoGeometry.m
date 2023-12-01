function [geo,angles]=readDiondoGeometry(folder)

% Written by P. Basford
% Based on Nikon data loader developed by A. Biguri and W. Sun
% W. Sun edited on 06.10.2018 for 3D pro version 5.2.6809.15380 (2018)

file = dir([folder, '\*.xml']);
if ~isempty(file)
    [geo,angles]=loadXMLDiondoGeometry(file);
end
end

function [geo,angles,height]=loadXMLDiondoGeometry(file)
    
    s = readstruct([file.folder,'\',file.name]);
    
    if s.ScanParameter.HelixMode ~= 'false'
        error("Helix mode reading not implemented")
    end
    if s.ScanParameter.Laminography ~= 'false'
        error("Laminography mode reading not implemented")
    end

    geo.nDetector = [s.Recon.ProjectionDimX; s.Recon.ProjectionDimY];
    geo.dDetector = [s.Recon.ProjectionPixelSizeX; s.Recon.ProjectionPixelSizeY];
    geo.sDetector = geo.nDetector.*geo.dDetector;

    geo.nVoxel = [s.Recon.VolumeDimX; s.Recon.VolumeDimY; s.Recon.VolumeDimZ];
    geo.dVoxel = [s.Recon.VolumeVoxelSizeX; s.Recon.VolumeVoxelSizeY;s.Recon.VolumeVoxelSizeZ];
    geo.sVoxel = geo.nVoxel.*geo.dVoxel;

    geo.DSD = s.Geometrie.SourceDetectorDist;
    geo.DSO = s.Geometrie.SourceObjectDist;
    geo.COR = s.Recon.ProjectionCenterOffsetX;

    if s.Recon.ProjectionCount ~= s.Recon.ProjectionCountPer360deg
        warning("Never tested with angles not in 360deg. If it doesnt work please contact tigre.toolbox@gmail.com to help improve this loader");
    end
    angles=linspace(0,2*pi,s.Recon.ProjectionCount);
    angles=-angles(1:end-1);


    height = char(s.Geometrie.Slices.string());
    height  = height(1:strfind(s.Geometrie.Slices.string,"mm")-2);
end
