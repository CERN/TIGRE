function [geo,angles]=from_vector_to_TIGRE(geo,vectors,varargin)
%FROM_VECTOR_TO_TIGRE TIGRE cone-beam geometry from per-view source and detector vectors.
%
%   [GEO, ANGLES] = FROM_VECTOR_TO_TIGRE(GEO, VECTORS) takes an N x 12 array
%   in the layout of ASTRA's cone_vec geometry (see FROM_TIGRE_TO_VECTOR)
%   and returns GEO with per-view DSO, DSD, offDetector, rotDetector,
%   offOrigin (zero) and COR (zero), plus dDetector and sDetector from the
%   lengths of u and v. ANGLES is 3 x N; pass it to Ax, Atb and the
%   algorithms. The image fields and nDetector are taken from the input GEO.
%
%   TIGRE's description of a view is redundant, so this returns one of many
%   equivalent geometries. Vector geometries no TIGRE geometry can describe
%   are refused with an error: u and v not perpendicular, a pixel size that
%   changes between views, a source at the centre of the image or in the
%   plane of the detector.
%
%   Options (name-value pairs):
%     'detector_reference'  'center' (default) or 'first_pixel': whether
%                           columns 4-6 hold the centre of the detector or
%                           of the first pixel.
%     'rtol'                Relative tolerance for the perpendicularity and
%                           pixel-size checks. Default 1e-6.
%
%   See also FROM_TIGRE_TO_VECTOR, AX, CHECKGEO.

[reference,rtol]=parseOptions(varargin{:});
if ~isnumeric(vectors) || ~ismatrix(vectors) || size(vectors,2)~=12 || isempty(vectors)
    error('TIGRE:from_vector_to_TIGRE:InvalidInput','vectors should be of size Nx12');
end
vectors=double(vectors);
if ~all(isfinite(vectors(:)))
    error('TIGRE:from_vector_to_TIGRE:InvalidInput','vectors must contain only finite values');
end
n=size(vectors,1);
nDetector=double(geo.nDetector(:));
scale=max(max(abs(vectors(:,1:6))));

DSO=zeros(1,n); DSD=zeros(1,n);
offDetector=zeros(2,n); rotDetector=zeros(3,n); angles=zeros(3,n);
dU=norm(vectors(1,7:9)); dV=norm(vectors(1,10:12));
if dU==0 || dV==0
    error('TIGRE:from_vector_to_TIGRE:InvalidInput','u and v must be non-zero in every view');
end
for i=1:n
    source=vectors(i,1:3).'; centre=vectors(i,4:6).'; u=vectors(i,7:9).'; v=vectors(i,10:12).';
    if strcmp(reference,'first_pixel')
        centre=centre+u*(nDetector(1)-1)/2+v*(nDetector(2)-1)/2;
    end
    if abs(norm(u)-dU)>rtol*dU || abs(norm(v)-dV)>rtol*dV
        error('TIGRE:from_vector_to_TIGRE:InvalidInput',...
            ['The length of u or v (the pixel size) changes between views (view %d). ',...
             'TIGRE uses one detector pixel size for the whole scan.'],i);
    end
    uHat=u/norm(u); vHat=v/norm(v);
    if abs(uHat.'*vHat)>rtol
        error('TIGRE:from_vector_to_TIGRE:InvalidInput',...
            'u and v are not perpendicular in view %d (angle %.6f degrees). TIGRE detectors have perpendicular axes.',...
            i,acosd(max(-1,min(1,uHat.'*vHat))));
    end
    vHat=vHat-(uHat.'*vHat)*uHat; vHat=vHat/norm(vHat);

    DSO(i)=norm(source);
    if DSO(i)<=rtol*scale
        error('TIGRE:from_vector_to_TIGRE:InvalidInput','The source must not be at the centre of the image (view %d)',i);
    end
    normal=cross(uHat,vHat);
    if abs((source-centre).'*normal)<=rtol*DSO(i)
        error('TIGRE:from_vector_to_TIGRE:InvalidInput','The source must not be in the plane of the detector (view %d)',i);
    end

    % Frame of the TIGRE angles: x axis towards the source; z axis as close
    % to world z as possible, so a source in the x-y plane gives (alpha,0,0).
    x=source/DSO(i);
    candidates=[[0;0;1] vHat uHat];
    for k=1:3
        z=candidates(:,k)-(candidates(:,k).'*x)*x;
        if norm(z)>1e-3
            z=z/norm(z);
            break
        end
    end
    angles(:,i)=eulerZYZ([x cross(z,x) z]);
    R=rz(angles(1,i))*ry(angles(2,i))*rz(angles(3,i));

    local=R.'*centre;
    DSD(i)=DSO(i)-local(1);
    offDetector(:,i)=[local(2); local(3)];
    Rd=R.'*[normal uHat vHat];
    [roll,pitch,yaw]=rollPitchYaw(Rd);
    rotDetector(:,i)=[yaw; pitch; roll];
end

geo.mode='cone';
geo.DSO=DSO;
geo.DSD=DSD;
geo.offDetector=offDetector;
geo.rotDetector=rotDetector;
geo.offOrigin=zeros(3,n);
geo.COR=zeros(1,n);
geo.dDetector=[dU; dV];
geo.sDetector=geo.nDetector(:).*geo.dDetector;
end

function [reference,rtol]=parseOptions(varargin)
reference='center'; rtol=1e-6;
if mod(numel(varargin),2)
    error('TIGRE:from_vector_to_TIGRE:InvalidInput','Options should be name-value pairs');
end
for k=1:2:numel(varargin)
    switch lower(varargin{k})
        case 'detector_reference'
            reference=lower(varargin{k+1});
            if ~any(strcmp(reference,{'center','first_pixel'}))
                error('TIGRE:from_vector_to_TIGRE:InvalidInput','detector_reference must be ''center'' or ''first_pixel''');
            end
        case 'rtol'
            rtol=varargin{k+1};
        otherwise
            error('TIGRE:from_vector_to_TIGRE:InvalidInput','Unknown option ''%s''',varargin{k});
    end
end
end

function a=eulerZYZ(R)
% Angles [alpha;theta;psi] with rz(alpha)*ry(theta)*rz(psi)==R.
theta=acos(max(-1,min(1,R(3,3))));
if sin(theta)>1e-9
    a=[atan2(R(2,3),R(1,3)); theta; atan2(R(3,2),-R(3,1))];
elseif R(3,3)>0
    a=[atan2(R(2,1),R(1,1)); theta; 0];
else
    a=[atan2(-R(2,1),R(2,2)); theta; 0];
end
end

function [roll,pitch,yaw]=rollPitchYaw(R)
% Angles with rz(roll)*ry(pitch)*rx(yaw)==R.
pitch=-asin(max(-1,min(1,R(3,1))));
if cos(pitch)>1e-9
    yaw=atan2(R(3,2),R(3,3));
    roll=atan2(R(2,1),R(1,1));
else
    yaw=0;
    roll=atan2(-R(1,2),R(2,2));
end
end

function R=rz(a)
R=[cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
end

function R=ry(a)
R=[cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)];
end
