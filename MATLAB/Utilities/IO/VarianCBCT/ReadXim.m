function img=ReadXim(filename, varargin)
%   READXIM(FILENAME, {READ_PIXEL_DATA})
%       reads XIM image files from Varian TrueBeam accelerators 
%       and returns a struct containing the XIM image information
%       and pixel data.
%   
%   FILENAME is the full path to the XIM-image file
%   READ_PIXEL_DATA determines whether pixel data should be decoded
%       0: skip
%       1: decode  (default)
%
%   Developed by Fredrik Nordström 2015
%   The program has not been tested with uncompressed images

read_pixel_data=1;
if nargin==2 && varargin{1}==0
    read_pixel_data=0;
end
fid = fopen(filename);

% Decode header
img.file_name=filename;
img.file_format_identifier=fread(fid,8,'*char')';
img.file_format_version=fread(fid,1,'*int32');
img.image_width=fread(fid,1,'*int32');
img.image_height=fread(fid,1,'*int32');
img.bits_per_pixel=fread(fid,1,'*int32');
img.bytes_per_pixel=fread(fid,1,'*int32');
img.compression_indicator=fread(fid,1,'*int32');

% Decode pixel data
if img.compression_indicator==1
    lookup_table_size = fread(fid,1,'int32');    
    lookup_table=fread(fid,lookup_table_size*4,'ubit2=>uint8');
    compressed_pixel_buffer_size = fread(fid,1,'int32');
    if read_pixel_data==1
        % Decompress image
        pixel_data=int32(zeros(img.image_width*img.image_height,1));
        pixel_data(1:img.image_width+1)=fread(fid,img.image_width+1,'*int32');
        lookup_table_pos=1;
        for image_pos=(img.image_width+2):(img.image_width*img.image_height)      
            if lookup_table(lookup_table_pos)==0
                diff=int32(fread(fid,1,'*int8'));
            elseif lookup_table(lookup_table_pos)==1
                diff=int32(fread(fid,1,'*int16'));
            else
                diff=int32(fread(fid,1,'*int32'));
            end
            pixel_data(image_pos)=diff+pixel_data(image_pos-1)+pixel_data(image_pos-img.image_width)-pixel_data(image_pos-img.image_width-1);
            lookup_table_pos=lookup_table_pos+1;
        end
        if img.bytes_per_pixel==2
            img.pixel_data=int16(reshape(pixel_data,img.image_width,img.image_height))';
        else
            img.pixel_data=reshape(pixel_data,img.image_width,img.image_height)';
        end
    else
        fseek(fid,compressed_pixel_buffer_size,'cof');
    end
    uncompressed_pixel_buffer_size = fread(fid,1,'*int32');
else
    uncompressed_pixel_buffer_size = fread(fid,1,'*int32');
    if read_pixel_data==1        
        switch img.bytes_per_pixel
            case 1
                pixel_data=fread(fid,uncompressed_pixel_buffer_size,'*int8');
            case 2
                pixel_data=fread(fid,uncompressed_pixel_buffer_size/2,'*int16');
            otherwise
                pixel_data=fread(fid,uncompressed_pixel_buffer_size/4,'*int32');
        end
        img.pixel_data=reshape(pixel_data,img.image_width,img.image_height)';
    else
        fseek(fid,uncompressed_pixel_buffer_size,'cof');
    end
end

% Decode histogram
number_of_bins_in_histogram = fread(fid,1,'*int32');
if number_of_bins_in_histogram>0
    img.histogram.number_of_bins_in_histogram=number_of_bins_in_histogram;
    img.histogram.histogram_data = fread(fid,number_of_bins_in_histogram,'*int32');
end

% Decode properties
number_of_properties = fread(fid,1,'*int32');
if number_of_properties>0
    img.properties=[];
end
for property_nr=1:number_of_properties
    property_name_length = fread(fid,1,'*int32');
    property_name = fread(fid,property_name_length,'*char')';
    property_type = fread(fid,1,'*int32');
    switch property_type
        case 0
            property_value = fread(fid,1,'*int32');
        case 1
            property_value = fread(fid,1,'double');
        case 2
            property_value_length = fread(fid,1,'*int32');
            property_value = fread(fid,property_value_length,'*char')';
        case 4
            property_value_length = fread(fid,1,'*int32');
            property_value = fread(fid,property_value_length/8,'double');
        case 5
            property_value_length = fread(fid,1,'*int32');
            property_value = fread(fid,property_value_length/4,'*int32');
        otherwise
            disp(' ')
            disp([property_name ': Property type ' num2str(property_type) ' is not supported! Aborting property decoding!']);
            fclose(fid);
            return;
    end
    img.properties=setfield(img.properties,property_name,property_value);
end

fclose(fid);