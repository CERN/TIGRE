function  showGeoCBCTDiagram()
%SHOWGEODIAGRAM Shows an image describing the Geometry of TIGRE

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

