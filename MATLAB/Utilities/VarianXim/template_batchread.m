
Ximfolder = 'E:\BigData\Edge\CBCT_Export\2019-09-03_121026';

filestr = dir([Ximfolder filesep '**' filesep 'Proj_*.xim']);

proj = [];
angle = [];
for ii = 1:length(filestr)
	filename = fullfile(filestr(ii).folder, filestr(ii).name);
	[page, rtn] = mexReadXim(filename);
	if(~isempty(page))
		proj(:,:,ii) = page;
		angle(ii) = rtn;
        ximstr{ii} = ReadXim(filename, 0);
	end
end
