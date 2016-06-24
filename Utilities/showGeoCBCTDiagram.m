function  showGeoCBCTDiagram()
%SHOWGEODIAGRAM Shows an image describing the Geometry of TIGRE
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
% 
% Copyright (c) 2015, University of Bath and 
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License:            Open Source under BSD. 
%                     See the full license at
%                     https://github.com/CERN/TIGRE/license.txt
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri
%--------------------------------------------------------------------------
if haveInet
figure('Name','Diagram of TIGRE Geometry');
title('Diagram of TIGRE Geometry');
geoimg=imread('https://i.imgur.com/mRweux3.png');
imshow(geoimg);

h = xlabel(''); 
pos = get(h,'Position'); 
delete(h)
h = title(char('Geometry definition for CBCT','    ©TIGRE toolbox','   DOI: XXX-XXXX XXXX'));
set(h,'Position',pos);
set(gca, 'XAxisLocation','top')
set(gcf, 'Color','white')

else
    disp('showGeoCBCTDiagram() needs Internet to work. Run doc(''TIGRE/Geometry'') to see the diagram.')
    
end

