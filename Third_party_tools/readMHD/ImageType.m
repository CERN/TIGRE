classdef ImageType < handle
%IMAGETYPE  Base class to handle nD images
%    This class provides basic functionality and structure for medical images, inspired by the Insight Toolkit (www.itk.org). It provides means to handle origin, spacing, orientation, etc.
%
%    Most of the functionality is adapted to n-D images, and has been tested for 2D, 3D and 4D images.
%
%   Example:
%
%     img = ImageType(size,origin,spacing,orientation);
%     img.data = rand(size);
%     img.GetPosition([1 2 1]'); % retrieve the world coordinates of voxel with index [1 2 1]'
%
%   See also ECHOIMAGETYPE, VECTORIMAGETYPE, MESHTYPE.

%   Written by Alberto Gomez 2011
%   King's College London
%   OpenSource code under the BSD-2 license
%   This software is distributed with no warranties.
%   This software was designed for research purposes and not for clinical use.
    
    properties(GetAccess = 'public', SetAccess = 'public')
        data=[];
        size=[];
        origin=[]; % first voxel is indexed by [1 1 1] !!!
        spacing= [];
        orientation = []; % orientation matrix (3 x 3)
        D = [];
        paddingValue = 0;
        MAX_CHUNK_SIZE = 150; % for internal operations, to preserve memory, maximum operation block will be of 100 x 100 x 100
    end
    
    methods(Access = public)
        %constructor: providing a standard constructor and a copy constructor
        function obj = ImageType(size,origin,spacing,orientation)
            if (nargin==1)
                % copy constructor. 
 
                obj = ImageType(size.size,size.origin,size.spacing,size.orientation);                 % argument "size" is another images
                obj.data = size.data;
                
            elseif(nargin > 0)
                obj.size = size(:);
                obj.spacing =  spacing(:);
                obj.origin = origin(:);
                obj.orientation = orientation;
                s = obj.size;
                obj.data = zeros(s(:)');
                
                obj.D = orientation;
                
                for i=1:numel(obj.size)
                    obj.D(:,i)=obj.D(:,i)*obj.spacing(i);
                end
                
            end
            
            
        end
        
        function [P, I] = GetBounds(obj,th,addBorder)
            % get the limits in x, y and z, in world coordinates (P) and in voxel coordinates (I)
            % Arguments 'th' , and 'addBorder' are optional. If th is specified, data below that threshold is ignored when calculating the bounds.
            % addBorder is a boolean which allows to add a 1 voxel border around the bounds.
            if (nargin==1)
                P = zeros(2*ndims(obj.data),1);

                
                c = zeros(ndims(obj.data),2^ndims(obj.data));
                
                  str1 = '';
                  str2 = '';
                  str3 = '';
                  intervals=[];
                  for i = 1:ndims(obj.data)
                        str1 = [str1 '[1 obj.size('  num2str(i) ')],'];
                        str2 = [str2 'i' num2str(i) ' '];
                        str3 = [str3 'i' num2str(i) '(:) '];
                  end
                  str1=str1(1:end-1); % remove last comma
                  eval([ '[' str2 ']= ndgrid(' str1 '); intervals = [' str3 ']; clear ' str2 ';']);
                
                for i = 1:2^ndims(obj.data)
                    c(:,i) = obj.GetPosition(intervals(i,:)');
                end
                
                P= reshape([min(c,[],2) max(c,[],2)]',1,[]);
                
               
                I = [-1 -1 -1 -1 -1 -1];
            elseif (nargin>1)
                % get the limits in x, y and z
                % get limit in x, y, z where abs value of intensity is over
                % the threshold th
                % project into xy;
                for i=1:ndims(obj.data)
                    otherdims = setdiff(1:ndims(obj.data),i);
                    data = abs(obj.data);
                    for j=otherdims
                        data = max(data,[],j);
                    end
                    bounds(i,1)=min(find(data));
                    bounds(i,2)=max(find(data));
                end
                
               
                P = obj.GetPosition(bounds);
                order = reshape([1:ndims(obj.data); (1:ndims(obj.data))+ndims(obj.data)],2*ndims(obj.data),[]);
                P=P(order);
                I = bounds(order);
                if (nargin==3)
                    if (addBorder)
                        I = I + reshape([-ones(1,ndims(obj.data)) ; ones(1,ndims(obj.data))],2*ndims(obj.data),[]); % leave a border
                        P = P + reshape([-obj.spacing' ; obj.spacing'],2*ndims(obj.data),[]); % leave a border
                    end
                end
            end
            P = P(:);
            I=I(:);
        end
        
        
        function pos = GetPosition(obj, index)
            % If obj is a K x n image, where K is the dimension (eg 2, 3 4, ...)
            % get the position of a voxel at index "index", a K x n matrix
            % also returns a K x n matrix
            
            nindices = size(index,2);
            pos = zeros(size(index));
            if nindices > obj.MAX_CHUNK_SIZE^size(index,1)
                % Divide into chunks
                NCHUNKS = ceil(nindices/( obj.MAX_CHUNK_SIZE^size(index,1)));
                chunked_size = ceil(nindices/NCHUNKS);
                intervals= 0:NCHUNKS-1;
                
                for i=1:numel(intervals)
                    ranges(1) = intervals(i)*chunked_size+1;
                    ranges(2) = min([(intervals(i)+1)*chunked_size ; nindices]);
                    nranges = ranges(2)-ranges(1)+1;
                    pos(:,ranges(1):ranges(2)) = obj.origin(:) * ones(1,nranges) + obj.D *(index(:,ranges(1):ranges(2))-1);
                end
                
            else
                pos = obj.origin(:) * ones(1,nindices) + obj.D *(index-1);
            end
        end
        
        function index = GetContinuousIndex(obj, pos)
            % get the continuous index position
            % pos has to be a K x n matrix where K is the dimensionality of
            % the image
            
            index =  obj.D\(pos - obj.origin(:)*ones(1,size(pos,2)) ) + 1;
        end
        
        function P = GetPixel(obj, index)
            % get the pixel value at an index position. "index" is a 3 x 1 matrix
            P = obj.data(index(1),index(2),index(3));
        end
        
        function c_in = isOutOfRange(obj,index)
            %  P = isInRange(index)
            % index must be a 3XN array
            % returns 0 if the index is inside the image, 0 otherwise
            
            
            c_in = zeros(1,size(index,2));
            ndims = numel(obj.size);
            for i=1:ndims
                c_in = c_in | index(i,:)<1 | index(i,:)>obj.size(i);
            end
            
        end
        
        function P = GetValue(obj, pos, varargin)
            % get the pixel value at a non index position
            % can use different interpolation schemes:
            % im.GetValue(pos)  returns the value using nearest neighrbor
            % interpolation, with pos given in world coordinates
            % im.GetValue(pos, mode) uses the following interpolation
            %   mode = 'NN'     nearest neighbor
            %   mode = 'linear'    (tri) linear interpolation
            %   mode = 'spline'    (cubic) b-spline interpolation
            
            P=obj.paddingValue;
            mode = 'NN'; % NN
            ndims = numel(obj.size);
            if  (size(varargin,2)>0)
                mode = varargin{1};
            end
            
            
            
            index = obj.GetContinuousIndex(pos);
            
            
            round_index = round(index);
            
            % find the indexes inside the range
            
            c_in = obj.isOutOfRange(round_index);
            
            if numel(c_in~=0)
                
                
                if (strcmp(mode,'NN'))
                    round_index(:,c_in) = round_index(:,c_in)*0+1;
                    
                    str = '';
                    
                    for i = 1:ndims
                        str = [str ', round_index('  num2str(i) ',:)'];
                    end
                    
                    in_1D = eval(['sub2ind(obj.size''' str ')']);
                    P = obj.data(in_1D);
                    
                elseif (strcmp(mode,'linear'))
                    
                    index(:,c_in) = index(:,c_in)*0+1;
                    P = obj.evaluate_at_point_linear(index);
                elseif (strcmp(mode,'spline'))
                    index(:,c_in) = index(:,c_in)*0+1;
                    P = obj.evaluate_at_point_spline(index);
                end
                P(c_in)=0;
                
            end
            
            
            
            
        end
        
        function SetOriginToCenter(obj)
            
            obj.origin = obj.orientation * ((obj.size -1).*obj.spacing/2.0 );
        end
        
        function out = extractFrame(obj,nframe)
            if numel(obj.size)~=4
                disp('WARNING: input image is not 4D. There might be problems')
            end
            out = ImageType(obj.size(1:end-1),obj.origin(1:end-1),obj.spacing(1:end-1),obj.orientation(1:end-1,1:end-1));
            out.data = obj.data(:,:,:,nframe);
            
        end
        
    end
    
    methods(Access = private)
        function  value = evaluate_at_point_linear(obj,continuous_index)
            % continuous_index is 3 x N
            
            %ind = find(obj.size<2);
            
            
            str = '';
            for i = 1:ndims(obj.data)
                str = [str ', continuous_index('  num2str(i) ',:)'];
            end
            
            value = eval(['interpn(double(obj.data)' str ',''linear'')']);
        end
        
        function  value = evaluate_at_point_spline(obj,continuous_index)
            % continuous_index is K x N
            str = '';
            for i = 1:ndims(obj.data)
                str = [str ', continuous_index('  num2str(i) ',:)'];
            end
            
            value = eval(['interpn(double(obj.data)' str ',''cubic'')']);
        end
        
    end % private methods
end
