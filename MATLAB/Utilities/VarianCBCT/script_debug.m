
datafolder = 'E:\BigData\Edge\CBCT_Export\2019-09-03_115522';

filestr = dir([datafolder filesep '**' filesep 'Proj_*.xim']);
%% load proj, XML
[proj, angle, blk, projinfo] = BatchReadXim(datafolder);
[geo, ScanXML, ReconXML] = GeometryFromXML(datafolder);

%% blk info
blkinfo = ReadXim('FilterBowtie.xim',0);

%% 
for ii = 1:length(angle)
    kvsrcrtn(ii) = projinfo{ii}.properties.KVSourceRtn;
    grtn(ii) = projinfo{ii}.properties.GantryRtn;
end

%%
slice = 80;
FDK_bowtie(FDK_bowtie<0) =0;
figure(1), imshow(FDK_bowtie(:,:,slice)/0.1)
figure(2), imshow(FDK_free(:,:,slice)/0.1)

