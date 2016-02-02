function [img info]=read_mhd(filename)
% This function is based upon "read_mhd" function from the package
% ReadData3D_version1 from the matlab exchange.
% Copyright (c) 2010, Dirk-Jan Kroon
% [image info ] = read_mhd(filename)

info = mhareadheader(filename);

[path name extension] = fileparts(filename);

if (isfield(info,'ElementNumberOfChannels'))
    ndims = str2num(info.ElementNumberOfChannels);
else
    ndims = 1;
end
img = ImageType();


if (ndims == 1)
    data = read_raw([ path  filesep  info.DataFile ], info.Dimensions,info.DataType,'native',0,ndims);
    img=ImageType(size(data),info.Offset',info.PixelDimensions',reshape(info.TransformMatrix,numel(info.PixelDimensions),numel(info.PixelDimensions)));
    img.data = data;
elseif (ndims == 3)
    clear img;
        
    [datax datay dataz] = read_raw([ path  filesep  info.DataFile], info.Dimensions,info.DataType,'native',0,ndims);
    img = VectorImageType(size(datax),info.Offset',info.PixelDimensions',reshape(info.TransformMatrix,numel(info.PixelDimensions),numel(info.PixelDimensions)));
    
    img.datax = datax; clear datax;
    img.datay=datay; clear datay;
    img.dataz = dataz; clear dataz;
    img.data = img.datax.^2+img.datay.^2+img.dataz.^2;
end

 

end

function [rawData rdy rdz] =read_raw(filename,imSize,type,endian,skip,ndims)

% Reads a raw file
% Inputs: filename, image size, image type, byte ordering, skip
% If you are using an offset to access another slice of the volume image
% be sure to multiply your skip value by the number of bytes of the 
% type (ie. float is 4 bytes).
% Inputs: filename, image size, pixel type, endian, number of values to
% skip.
% Output: image

rdy=[];
rdz=[];
fid = fopen(filename,'rb',endian);
if (fid < 0)
    display(['Filename ' filename ' does not exist']);
    rawData = -1;
else
    if (ndims == 1)
        status = fseek(fid,skip,'bof');
        if status == 0
            rawData = fread(fid,prod(imSize),type);
            fclose(fid);
            rawData = reshape(rawData,imSize);
        else
            rawData = status;
        end
    else
        %disp('Vector Image');
        status = fseek(fid,skip,'bof');
        if status == 0
            r = fread(fid,prod(imSize)*3,type);
            fclose(fid);
            if length(imSize) == 3
            %    slices = length(r)/imSize(1)/imSize(2)/3;% 3 for vx, vy, vz
            %    imSize(3) = slices;
            %    imSize(4) = 3; 
                im_size=[ 3 imSize([1 2 3]) ];
            elseif length(imSize) == 4
                im_size=[3 imSize([1 2 3 4]) ];
            end
            r = reshape(r,im_size);
           
            
        else
            r = status;
        end
         if length(imSize) == 3
            rawData=squeeze(r(1,:,:,:));
            rdy=squeeze(r(2,:,:,:));
            rdz=squeeze(r(3,:,:,:));
         elseif length(imSize) == 4
             rawData=squeeze(r(1,:,:,:,:));
            rdy=squeeze(r(2,:,:,:,:));
            rdz=squeeze(r(3,:,:,:,:));
         end
        

    end
end

end



function info =mhareadheader(filename)
% Function for reading the header of a Insight Meta-Image (.mha,.mhd) file
% 
% info  = mha_read_header(filename);
%
% examples:
% 1,  info=mha_read_header()
% 2,  info=mha_read_header('volume.mha');
info = [];
if(exist('filename','var')==0)
    [filename, pathname] = uigetfile('*.mha', 'Read mha-file');
    filename = [pathname filename];
end

fid=fopen(filename,'rb');
if(fid<0)
    fprintf('could not open file %s\n',filename);
    return
end

info.Filename=filename;
info.Format='MHA';
info.CompressedData='false';
readelementdatafile=false;
while(~readelementdatafile)
    str=fgetl(fid);
    s=find(str=='=',1,'first');
    if(~isempty(s))
        type=str(1:s-1); 
        data=str(s+1:end);
        while(type(end)==' '); type=type(1:end-1); end
        while(data(1)==' '); data=data(2:end); end
    else
        type=''; data=str;
    end
    
    switch(lower(type))
        case 'ndims'
            info.NumberOfDimensions=sscanf(data, '%d')';
        case 'dimsize'
            info.Dimensions=sscanf(data, '%d')';
        case 'elementspacing'
            info.PixelDimensions=sscanf(data, '%lf')';
        case 'elementsize'
            info.ElementSize=sscanf(data, '%lf')';
            if(~isfield(info,'PixelDimensions'))
                info.PixelDimensions=info.ElementSize;
            end
        case 'elementbyteordermsb'
            info.ByteOrder=lower(data);
        case 'anatomicalorientation'
            info.AnatomicalOrientation=data;
        case 'centerofrotation'
            info.CenterOfRotation=sscanf(data, '%lf')';
        case 'offset'
            info.Offset=sscanf(data, '%lf')';
        case 'binarydata'
            info.BinaryData=lower(data);
        case 'compresseddatasize'
            info.CompressedDataSize=sscanf(data, '%d')';
        case 'objecttype',
            info.ObjectType=lower(data);
        case 'transformmatrix'
            info.TransformMatrix=sscanf(data, '%lf')';
        case 'compresseddata';
            info.CompressedData=lower(data);
        case 'binarydatabyteordermsb'
            info.ByteOrder=lower(data);
        case 'elementdatafile'
            info.DataFile=data;
            readelementdatafile=true;
        case 'elementtype'
            info.DataType=lower(data(5:end));
        case 'headersize'
            val=sscanf(data, '%d')';
            if(val(1)>0), info.HeaderSize=val(1); end
        otherwise
            info.(type)=data;
    end
end

switch(info.DataType)
    case 'char', info.BitDepth=8;
    case 'uchar', info.BitDepth=8;
    case 'short', info.BitDepth=16;
    case 'ushort', info.BitDepth=16;
    case 'int', info.BitDepth=32;
    case 'uint', info.BitDepth=32;
    case 'float', info.BitDepth=32;
    case 'double', info.BitDepth=64;
    otherwise, info.BitDepth=0;
end
if(~isfield(info,'HeaderSize'))
    info.HeaderSize=ftell(fid);
end
fclose(fid);
end