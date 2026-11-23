function vectors=from_TIGRE_to_vector(geo,angles,varargin)
%FROM_TIGRE_TO_VECTOR Per-view source and detector vectors of a TIGRE cone-beam geometry.
%
%   VECTORS = FROM_TIGRE_TO_VECTOR(GEO, ANGLES) returns an N x 12 array, one
%   row per view, in the layout of ASTRA's cone_vec geometry:
%
%       [srcX srcY srcZ  dX dY dZ  uX uY uZ  vX vY vZ]
%
%   src is the source position and d the centre of the detector. u is the
%   vector from one detector pixel to the next along U (the first entry of
%   geo.nDetector), and v the vector to the next pixel along V. The vectors
%   are in TIGRE's world frame, in the units of geo (normally mm), with the
%   centre of the image at the origin.
%
%   ANGLES is 1 x N or 3 x N, as for Ax. DSO, DSD, offOrigin, offDetector,
%   rotDetector and COR may be given once or per view. Only cone-beam
%   geometries are supported.
%
%   FROM_TIGRE_TO_VECTOR(..., 'detector_reference', 'first_pixel') puts the
%   centre of the first pixel, instead of the centre of the detector, in
%   columns 4-6.
%
%   See also FROM_VECTOR_TO_TIGRE, AX, CHECKGEO.

reference=parseReference(varargin{:});
if isfield(geo,'mode') && ~strcmp(geo.mode,'cone')
    error('TIGRE:from_TIGRE_to_vector:InvalidInput','Only cone-beam geometries are supported, got mode ''%s''',geo.mode);
end
angles=double(angles);
if size(angles,1)==1
    angles=[angles; zeros(2,size(angles,2))];
end
if size(angles,1)~=3 || isempty(angles)
    error('TIGRE:from_TIGRE_to_vector:InvalidInput','Angles should be of size 1xN or 3xN');
end
n=size(angles,2);

DSO=perView(geo,'DSO',1,n,[]);
DSD=perView(geo,'DSD',1,n,[]);
offOrigin=perView(geo,'offOrigin',3,n,zeros(3,1));      % (x;y;z)
offDetector=perView(geo,'offDetector',2,n,zeros(2,1));  % (u;v)
rotDetector=perView(geo,'rotDetector',3,n,zeros(3,1));  % about (x;y;z)
COR=perView(geo,'COR',1,n,0);
dDetector=double(geo.dDetector(:));                     % (u;v)
nDetector=double(geo.nDetector(:));

vectors=zeros(n,12);
for i=1:n
    R=rotZYZ(angles(1,i),angles(2,i),angles(3,i));
    Rd=rotRollPitchYaw(rotDetector(3,i),rotDetector(2,i),rotDetector(1,i));
    shift=COR(i)*[-sin(angles(1,i)); cos(angles(1,i)); 0]-offOrigin(:,i);
    source=R*[DSO(i);0;0]+shift;
    centre=R*[DSO(i)-DSD(i); offDetector(1,i); offDetector(2,i)]+shift;
    F=R*Rd;
    u=dDetector(1)*F(:,2);
    v=dDetector(2)*F(:,3);
    if strcmp(reference,'first_pixel')
        centre=centre-u*(nDetector(1)-1)/2-v*(nDetector(2)-1)/2;
    end
    vectors(i,:)=[source.' centre.' u.' v.'];
end
end

function reference=parseReference(varargin)
reference='center';
if isempty(varargin)
    return
end
if numel(varargin)~=2 || ~strcmpi(varargin{1},'detector_reference')
    error('TIGRE:from_TIGRE_to_vector:InvalidInput','The only option is ''detector_reference''');
end
reference=lower(varargin{2});
if ~any(strcmp(reference,{'center','first_pixel'}))
    error('TIGRE:from_TIGRE_to_vector:InvalidInput','detector_reference must be ''center'' or ''first_pixel''');
end
end

function value=perView(geo,name,rows,n,default)
if isfield(geo,name)
    value=double(geo.(name));
elseif ~isempty(default)
    value=default;
else
    error('TIGRE:from_TIGRE_to_vector:InvalidInput','geo.%s is not set',name);
end
if numel(value)==rows
    value=repmat(value(:),1,n);
elseif rows==1 && numel(value)==n
    value=value(:).';
elseif ~isequal(size(value),[rows n])
    error('TIGRE:from_TIGRE_to_vector:InvalidInput','geo.%s should be %dx1 or %dx%d',name,rows,rows,n);
end
end

function R=rotZYZ(alpha,theta,psi)
R=rz(alpha)*ry(theta)*rz(psi);
end

function R=rotRollPitchYaw(roll,pitch,yaw)
R=rz(roll)*ry(pitch)*rx(yaw);
end

function R=rz(a)
R=[cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
end

function R=ry(a)
R=[cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)];
end

function R=rx(a)
R=[1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
end
