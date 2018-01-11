function [geo,angles,whitelevel]=readXtekctGeometry(fpath)

mylist=cellstr(ls(fpath));
filexct=find(~cellfun('isempty', strfind(mylist,'.xtekct'))',1);
filexct=mylist{filexct};


fid=fopen([fpath,'\',filexct]);
if(fid==-1)
   error('Wrong file path');
end
xtekctText = textscan(fid, '%s %s', 'Delimiter', '=', 'HeaderLines', 1, 'CommentStyle', '[');
fclose(fid);

%% Detector information
% Number of pixel in the detector
geo.nDetector=[str2double(xtekctText{2}(strcmp('DetectorPixelsX', xtekctText{1})));
               str2double(xtekctText{2}(strcmp('DetectorPixelsY', xtekctText{1})))  ];
% Size of pixels in the detector       
geo.dDetector=[str2double(xtekctText{2}(strcmp('DetectorPixelSizeX', xtekctText{1})));
               str2double(xtekctText{2}(strcmp('DetectorPixelSizeY', xtekctText{1})))  ];
% Total size of the detector
geo.sDetector=geo.nDetector.*geo.dDetector;

%% Offset of the detector:
% NOTE: This is untested. I do not have Nikkon's Geometry definition, so I
% dont know if:
% 1-The direction (sign) of these matches with TIGRE's geometry definition
% 2-The values are in mm or pixels
%
% Please contact tigre.toolbox@gmail.com if this doesnt work/you have more
% information
geo.offDetector=[str2double(xtekctText{2}(strcmp('DetectorOffsetX', xtekctText{1})));
                 str2double(xtekctText{2}(strcmp('DetectorOffsetY', xtekctText{1})))  ];
           
%% Image information
% Number of pixel in the detector
geo.nVoxel=[str2double(xtekctText{2}(strcmp('VoxelsX', xtekctText{1})));
               str2double(xtekctText{2}(strcmp('VoxelsY', xtekctText{1})));
               str2double(xtekctText{2}(strcmp('VoxelsZ', xtekctText{1}))) ];
% Size of each pixel           
geo.dVoxel=[str2double(xtekctText{2}(strcmp('VoxelSizeX', xtekctText{1})));
               str2double(xtekctText{2}(strcmp('VoxelSizeY', xtekctText{1})));
               str2double(xtekctText{2}(strcmp('VoxelSizeZ', xtekctText{1}))) ];
% Size of the image in mm           
geo.sVoxel=geo.nVoxel.*geo.dVoxel;
geo.offOrigin=[0;0;0];
warning('The image size has been untouched. Often this is bigger than what TIGRE v1.1.X can handle. Consider changing the number of voxels.')
%% Global geometry
geo.DSO=str2double(xtekctText{2}(strcmp('SrcToObject', xtekctText{1})));
geo.DSD=str2double(xtekctText{2}(strcmp('SrcToDetector', xtekctText{1})));
geo.COR=-str2double(xtekctText{2}(strcmp('CentreOfRotationTop', xtekctText{1})));
if (geo.COR==0)
    warning('Centre of Rotation seems to be zero. Make sure that it is true and that the machine did not omit that information. Consider computing the COR with computeCOR() function.');
else
    warning('TIGRE doesnt know if the sign of COR is the right one. Consider triying both and reporting to tigre.toolbox@gmail.com.');
end
%% whitelevel

whitelevel=str2double(xtekctText{2}(strcmp('WhiteLevel', xtekctText{1})));

%% angles
filexct=find(~cellfun('isempty', strfind(mylist,'.ang'))',1);
filexct=mylist{filexct};
fid=fopen([fpath,'\',filexct]);
xtekctText = textscan(fid, '%s %s', 'Delimiter', '\t', 'HeaderLines', 1);
fclose(fid);
angles=str2double(xtekctText{2})*pi/180;
angles=angles';


end