clear,clc
%%
folder = 'E:\BigData\Edge\CBCT_Export\2020-03-10_103712\Acquisitions\789186315';
filestr = dir(fullfile(folder, '*.xim'));

tic
for ii = 1:length(filestr)
    ii
    filename = fullfile(folder, filestr(ii).name);
    XIM{ii} = ReadXim(filename, 0);
end
toc

%%
for ii=1:length(XIM)
    kvCurve(ii) = XIM{ii}.properties.KVKiloVolts;
    kvms(ii) = XIM{ii}.properties.KVMilliSeconds;
    kvmA(ii) = XIM{ii}.properties.KVMilliSeconds;
    kvPulseCounter(ii) = XIM{ii}.properties.KVBeamPulseCounter;
end

