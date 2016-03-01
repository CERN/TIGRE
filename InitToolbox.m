% Add tolbox folders
addpath('.\Algorithms');
addpath('.\Utilities');
addpath('.\Test_data');
addpath('.\Mex_files');
addpath('.\Utilities\Quality_measures');

% Perceptually uniform colormaps
addpath('.\Colormaps');
% Add third party tools from FEX
addpath('.\Third_party_tools\arrow3d'); % 3D shepp-Logan
addpath('.\Third_party_tools\sec2hours'); 
addpath('.\Third_party_tools\readMHD'); 


[user, sys]=memory;

if sys.PhysicalMemory.Total<9000000000 % 8Gb
   warning('Your Computer has 8Gb or less of RAM memory. Using image sizes of higher than 512^3 is not recomended') 
end

if sys.PhysicalMemory.Total<2500000000 % 8Gb
   warning('Your Computer has 2Gb or less of RAM memory. Using image sizes of higher than 256^3 is not recomended') 
end