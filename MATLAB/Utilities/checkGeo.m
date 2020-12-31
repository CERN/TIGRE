function geo=checkGeo(geo, angles )
% CHECKGEO checks for correctness in the geometry struct and prepares it
% for usage in the mex files. 
geofields_mandatory={'nVoxel','sVoxel','dVoxel', ...
    'nDetector','sDetector','dDetector',...
    'DSO','DSD'};

geofields_optional={'offOrigin','offDetector','rotDetector','COR',...
    'mode','accuracy'};
allfields=horzcat(geofields_mandatory,geofields_optional);
fnames=fieldnames(geo);
% Find if geo has fields we do not recongize
unknown=~ismember(fnames,allfields);
% there must be not unknown variables
% TODO: Do we really want to enforce this? Perhaps just a warning?
assert(~sum(unknown),'TIGRE:checkGeo:BadGeometry',['The following fields are not known by TIGRE:\n' strjoin(fnames(unknown)),'\nMake sure you have not misspelled any field or introduced unnecesary fields.'])



% check mandatory fields
mandatory=ismember(geofields_mandatory,fnames);
assert(sum(mandatory)==length(geofields_mandatory),'TIGRE:checkGeo:BadGeometry',['The following fields are missing:\n' strjoin(geofields_mandatory(~mandatory))])

% Lets check the validity of them now. We need to be explicit in there

%check if they are double
for ii=1:length(geofields_mandatory)
    assert(isa(geo.(geofields_mandatory{ii}),'double'),'TIGRE:checkGeo:BadGeometry',['Field geo.', geofields_mandatory{ii},' is not double type.'])
end

% Image data
assert(isequal(size(geo.nVoxel),[3 1]),'TIGRE:checkGeo:BadGeometry','geo.nVoxel should be 3x1')
assert(isequal(geo.nVoxel,round(geo.nVoxel)),'TIGRE:checkGeo:BadGeometry','geo.nVoxel should be a natural number.')

assert(isequal(size(geo.sVoxel),[3 1]),'TIGRE:checkGeo:BadGeometry','geo.sVoxel should be 3x1')

assert(isequal(size(geo.dVoxel),[3 1]),'TIGRE:checkGeo:BadGeometry','geo.sVoxel should be 3x1')

assert(sum(abs(geo.dVoxel.*geo.nVoxel-geo.sVoxel))<1e-6, 'TIGRE:checkGeo:BadGeometry', 'nVoxel*dVoxel is not sVoxel, something is wrong in the numbers')

% Detector data
assert(isequal(size(geo.nDetector),[2 1]),'TIGRE:checkGeo:BadGeometry','geo.nDetector should be 2x1')
assert(isequal(geo.nDetector,round(geo.nDetector)),'TIGRE:checkGeo:BadGeometry','geo.nDetector should be a natural number.')

assert(isequal(size(geo.sDetector),[2 1]),'TIGRE:checkGeo:BadGeometry','geo.sDetector should be 2x1')

assert(isequal(size(geo.dDetector),[2 1]),'TIGRE:checkGeo:BadGeometry','geo.sDetector should be 2x1')

assert(sum(abs(geo.dDetector.*geo.nDetector-geo.sDetector))<1e-6, 'TIGRE:checkGeo:BadGeometry', 'nDetector*dDetector is not sDetector, something is wrong in the numbers')

% DSD DSO

assert(isequal(size(geo.DSD),[1 1]) | isequal(size(geo.DSD),[1 size(angles,2)]),'TIGRE:checkGeo:BadGeometry','geo.DSD Should be 1x1 or 1xsize(angles,2)')
if isequal(size(geo.DSD),[1 1])
    geo.DSD=repmat(geo.DSD,[1, size(angles,2)]);
end
assert(isequal(size(geo.DSO),[1 1]) | isequal(size(geo.DSO),[1 size(angles,2)]),'TIGRE:checkGeo:BadGeometry','geo.DSD Should be 1x1 or 1xsize(angles,2)')

if isequal(size(geo.DSO),[1 1])
    geo.DSO=repmat(geo.DSO,[1, size(angles,2)]);
    
end

assert(all(geo.DSD>=geo.DSO), 'TIGRE:checkGeo:BadGeometry','DSD should be bigger or equal to DSO');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we know that optional fields are properly written or they would have
% been flagged before.
% We need to be explicit here again.
% geofields_optional={'offOrigin','offDetector','rotDetector','COR',...
%                     'mode','accuracy'};

if isfield(geo,'offOrigin')
    assert(isequal(size(geo.offOrigin),[3 1]) | isequal(size(geo.offOrigin),[3 size(angles,2)]),'TIGRE:checkGeo:BadGeometry','geo.offOrigin Should be 3x1 or 3xsize(angles,2)')
    assert(isa(geo.offOrigin,'double'),'TIGRE:checkGeo:BadGeometry','Field geo.offOrigin is not double type.' )
    if isequal(size(geo.offOrigin),[3 1])
        geo.offOrigin=repmat(geo.offOrigin,[1, size(angles,2)]);
    end
else
    geo.offOrigin=zeros(3,size(angles,2));
end

if isfield(geo,'offDetector')
    assert(isequal(size(geo.offDetector),[2 1]) | isequal(size(geo.offDetector),[2 size(angles,2)]),'TIGRE:checkGeo:BadGeometry','geo.offDetector Should be 2x1 or 2xsize(angles,2)')
    assert(isa(geo.offDetector,'double'),'TIGRE:checkGeo:BadGeometry','Field geo.offDetector is not double type.' )
    if isequal(size(geo.offDetector),[2 1])
        geo.offDetector=repmat(geo.offDetector,[1, size(angles,2)]);
    end
else
    geo.offDetector=zeros(2,size(angles,2));
end

if isfield(geo,'rotDetector')
    assert(isequal(size(geo.rotDetector),[3 1]) | isequal(size(geo.rotDetector),[3 size(angles,2)]),'TIGRE:checkGeo:BadGeometry','geo.rotDetector Should be 3x1 or 3xsize(angles,2)')
    assert(isa(geo.rotDetector,'double'),'TIGRE:checkGeo:BadGeometry','Field geo.rotDetector is not double type.' )
    if isequal(size(geo.rotDetector),[3 1])
        geo.rotDetector=repmat(geo.rotDetector,[1, size(angles,2)]);
    end
else
    geo.rotDetector=zeros(3,size(angles,2));
end

if isfield(geo,'COR')
    assert(isequal(size(geo.COR),[1 1]) | isequal(size(geo.COR),[1 size(angles,2)]),'TIGRE:checkGeo:BadGeometry','geo.COR Should be 1x1 or 1xsize(angles,2)')
    assert(isa(geo.COR,'double'),'TIGRE:checkGeo:BadGeometry','Field geo.COR is not double type.' )
    if isequal(size(geo.COR),[1 1])
        geo.COR=repmat(geo.COR,[1, size(angles,2)]);
    end
else
    geo.COR=zeros(1,size(angles,2));
end

if isfield(geo,'mode')
    assert(ischar(geo.mode),'TIGRE:checkGeo:BadGeometry','geo.mode should be a character array');
    assert(strcmp(geo.mode,'cone')|strcmp(geo.mode,'parallel'),'TIGRE:checkGeo:BadGeometry','geo.mode should ''cone'' or ''parallel''')
else
    geo.mode='cone';
end

if isfield(geo,'accuracy')
    assert(isscalar(geo.accuracy),'TIGRE:checkGeo:BadGeometry','geo.accuracy should be a scalar');
    assert(isa(geo.accuracy,'double'),'TIGRE:checkGeo:BadGeometry','geo.accuracy should be double');
    if geo.accuracy>1
        warning('geo.accuracy too big, you will ignore image information resulting in wrong reconstruction.\n Change geo.accuracy to smaller or equal than 1.')
    end
else
    geo.accuracy=0.5;
end


end

