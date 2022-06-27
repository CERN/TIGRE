function BHCalib = BH_SpectrumFilter(BHCalib)
%BH_SPECTRUMFILTER Summary of this function goes here
%   Detailed explanation goes here

% No filter at all
if(contains(BHCalib.filter.name,'None'))
    BHCalib.filter.spec = BHCalib.source.spec;
    return;
end
    
% filter thickness
thickness = BHCalib.filter.thickness;
% miu
ac = BHCalib.filter.ac;

% I/I_0 = exp(- length * miu): dependent on incident energy
atten = exp(-thickness.*ac);
BHCalib.filter.spec = BHCalib.source.spec .*atten(1:length(BHCalib.source.spec));

% normalized the filtered spectrum
BHCalib.filter.spec = BHCalib.filter.spec./max(BHCalib.filter.spec);

end

