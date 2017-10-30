classdef ImageType < handle
    % This class is compliant with itk::Image class
    % by Alberto Gomez, 2011
    %
    % Update: add Support for passing several points, for speed (start 18/oct/2011)
    
    properties(GetAccess = 'public', SetAccess = 'public')
        data=[];
        size=[];
        origin=[]; % first voxel is indexed by [1 1 1] !!!
        spacing= [];
        orientation = []; % orientation matrix (3 x 3)
        index = [];
        D = [];
        paddingValue = 0;
        MAX_CHUNK_SIZE = 200; % for internal operations, to preserve memory, maximum operation block will be of 100 x 100 x 100
        ndimensions = 3; % by default 3D
    end
    
    
    methods(Access = public)
        %constructor
        function obj = ImageType(size,origin,spacing,orientation)
            if (nargin==1)
                % argument "size" is another images
                obj = ImageType(size.size,size.origin,size.spacing,size.orientation);
                obj.index  =size.index;
            elseif (nargin==2) &&  strcmp(origin,'copy')
                obj = ImageType(size.size,size.origin,size.spacing,size.orientation);
                obj.data = size.data;
                obj.index  =size.index;
            elseif(nargin > 0)
                obj.size = size(:);
                obj.spacing =  spacing(:);
                obj.origin = origin(:);
                obj.orientation = orientation;
                s = obj.size;
                if numel(size)==1
                    obj.data = zeros(s(:)',1);
                    obj.size = size(:);
                else
                    obj.data = zeros(s(:)');
                end
                
                obj.D = orientation;
                
                for i=1:numel(obj.spacing)
                    obj.D(:,i)=obj.D(:,i)*obj.spacing(i);
                end
                obj.ndimensions = numel(obj.origin);
                obj.index = obj.spacing*0+1; % by default, the index points to the first pixel.
                
            end
            
            
        end
        
        function im = force2D(obj)
            im = ImageType(obj);
            im.data = obj.data;
            im.size= im.size(1:2);
            im.orientation = im.orientation(1:2,1:2);
            im.origin = im.origin(1:2);
            im.spacing= im.spacing(1:2);
            im.D = im.D(1:2,1:2);
        end
        
        function [P I] = GetBounds(obj,th,addBorder)
            % get the limits in x, y and z, in wc and in index
            %addborder is a boolean
            if (nargin==1)
                P = zeros(2*ndims(obj.data),1);
                % find the cornes
                nn = numel(obj.spacing);
                
                
                str1 = '';
                str2 = '';
                str3 = '';
                intervals=[];
                for i = 1:nn
                    str1 = [str1 '[1 obj.size('  num2str(i) ')],'];
                    str2 = [str2 'i' num2str(i) ' '];
                    str3 = [str3 'i' num2str(i) '(:) '];
                end
                str1=str1(1:end-1); % remove last comma
                eval([ '[' str2 ']= ndgrid(' str1 '); intervals = [' str3 ']; clear ' str2 ';']);
                intervals = unique(intervals,'rows')+((obj.index-1)*ones(1,size(unique(intervals,'rows'),1)))';
                %intervals = unique(intervals,'rows');%+((obj.index-1)*ones(1,size(unique(intervals,'rows'),1)))'
                c = zeros(nn,size(intervals,1));
                for i = 1:size(intervals,1)
                    c(:,i) = obj.GetPosition(intervals(i,:)');
                end
                
                P= reshape([min(c,[],2) max(c,[],2)]',1,[]);
                
                
                I = [-1 -1 -1 -1 -1 -1];
            elseif (nargin>1)
                % get the limits in x, y and z
                % get limit in x, y, z where abs value of intensity is over
                % the threshold th
                % project into xy;
                for i=1:numel(obj.size)
                    otherdims = setdiff(1:ndims(obj.data),i);
                    data = abs(obj.data);
                    for j=otherdims
                        data = max(data,[],j);
                    end
                    bounds(i,1)=min(find(data));
                    bounds(i,2)=max(find(data));
                end
                
                
                P = obj.GetPosition(bounds);
                order = reshape([1:numel(obj.size); (1:numel(obj.size))+numel(obj.size)],2*numel(obj.size),[]);
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
        
        % The index is the index starting at 1.
        function pos = GetPosition(obj, varargin)
            % If obj is a K x n image, where K is the dimension (eg 2, 3 4, ...)
            % get the position of a voxel at index "index", a K x n matrix
            % also returns a K x n matrix
            %
            
            add_im_index = false;
            if numel(varargin)==1
                index = varargin{1};
            else
                add_im_index = true;
                index = (1:prod(obj.size));
            end
            
            nindices = size(index,2);
            if size(index,1)==1 % If index is 1D, is assumes conversion from ind to sub
                index = obj.oned_to_nd_index(index);
%                 str = '';
%                 str2 = '';
%                 for i = 1:numel(obj.spacing)
%                     str = [str 'i'  num2str(i) ','];
%                     str2 = [str2 'i'  num2str(i) '(:) '];
%                 end
%                 
%                 eval([ '[' str(1:end-1) ']=ind2sub(obj.size'',index);']);
%                 % add the image index
%                 if add_im_index
%                     eval([ 'index =  ([' str2(1:end-1) ']''+(obj.index-1)*ones(1,numel(i1)));'  ]);
%                 else
%                     eval([ 'index =  [' str2(1:end-1) ']'';'  ]);
%                 end
            end
            
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
                    pos(:,ranges(1):ranges(2)) = obj.origin(:) * ones(1,nranges) + obj.D *(index(:,ranges(1):ranges(2))-ones(numel(obj.spacing),nranges));
                    %pos(:,ranges(1):ranges(2)) = obj.origin(:) * ones(1,nranges) + obj.D *(index(:,ranges(1):ranges(2))-1 +obj.index*ones(1,numel(ranges(1):ranges(2)))-1);
                end
                
            else
                pos = obj.origin(:) * ones(1,nindices) + obj.D *(index- ones(numel(obj.spacing),nindices)  );
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
            for i=1:obj.ndimensions
                c_in = c_in | (index(i,:)<obj.index(i)) | (index(i,:)>(obj.size(i)+obj.index(i)-1));
            end
            
        end
        
        function [N,s] = GetNeighbourhood(obj,varargin)
            % get the 8-connected neighbourhood (3D) or the K-connected
            % neighbourhood for other dimensions
            % s is a column vector with the extent of the resulting
            % neighbourhood
            
            neig = 1;
            if numel(varargin)
                neig = varargin{1};
            end
            
            if numel(neig) ==1
                neig = neig*ones(obj.ndimensions,1);
            end
            
            s = zeros(obj.ndimensions*2,1);
            str1='';
            str2='';
            str3='';
            for i=1:obj.ndimensions
                str1 = [str1 'x' num2str(i) ','];
                str2 = [str2 'x' num2str(i) '(:) '];
                str3 = [str3 '-' num2str(neig(i)) ':' num2str(neig(i)) ',' ];
                s((i-1)*2+[1 2]')=neig(i)*[-1 1];
            end
            eval([ '[' str1(1:end-1) ']=ndgrid(' str3(1:end-1) ');']); % [x,y,z]=ndgrid(-1:1,-1:1,-1:1);
            eval([ 'N = [' str2(1:end-1) '];' ]);
            
        end
        
        function N = GetNeighbourhoodBall(obj,varargin)
            % get the 4-connected neighbourhood (3D) or the K-connected
            % neighbourhood for other dimensions
            
            
            neig = 1;
            if numel(varargin)
                neig = varargin{1};
            end
            
            if numel(neig) ==1
                neig = neig*ones(obj.ndimensions,1);
            end
            
            
            str1='';
            str2='';
            str3='';
            str4 = '';
            for i=1:obj.ndimensions
                str1 = [str1 'x' num2str(i) ','];
                str2 = [str2 'x' num2str(i) '(:) '];
                str3 = [str3 '-' num2str(neig(i)) ':' num2str(neig(i)) ',' ];
                str4 = [str4 '(x' num2str(i) '/' num2str(neig(i)) ').^2+'];
                
            end
            
            
            eval([ '[' str1(1:end-1) ']=ndgrid(' str3(1:end-1) ');']); % [x,y,z]=ndgrid(-1:1,-1:1,-1:1);
            
            eval(['mask = ' str4(1:end-1) '<=1;']); % mask = (x/neigh(1)).^2 + (y/neigh(2)).^2 + (z/neigh(3)).^2 <= 1;
            eval([ 'N = [' str2(1:end-1) '];' ]);
            N(mask(:)==0,:)=[];
            
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
            
            P=ones(size(pos,2),1)*obj.paddingValue;
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
                    %round_index(:,c_in) = round_index(:,c_in)*0+1;
                    
                    in_1D=obj.nd_to_oned_index(round_index(:,~c_in));
                    P(~c_in) = obj.data(in_1D);
                    
                elseif (strcmp(mode,'linear'))
                    
                    index(:,c_in) = index(:,c_in)*0+1;
                    P = obj.evaluate_at_point_linear(index);
                elseif (strcmp(mode,'spline'))
                    index(:,c_in) = index(:,c_in)*0+1;
                    P = obj.evaluate_at_point_spline(index);
                end
                P(c_in)=obj.paddingValue;
                
            end
            
            
            
            
        end
        
       
        function ndindex = oned_to_nd_index(obj, oned_index_patch)
            
            total_size = obj.size';
            
            str = [];
            str2 =[];
            for ii=1:numel(obj.spacing)
                str = [str 'i' num2str(ii) ','];
                str2 = [str2 'i' num2str(ii) '(:) '];
            end
            eval(['['  str(1:end-1) '] = ind2sub(total_size,oned_index_patch); index = [' str2(1:end-1) ']'';' ]);
            ndindex = index+ obj.index*ones(1,size(index,2))-1;
        end
        
         % In this function, nd_index_patch must be, each index, in a column
        function onedindex = nd_to_oned_index(obj, nd_index_patch)
            
            if ~isprop(obj,'borderIndex1')
                index_local = nd_index_patch;
                total_size = obj.size';     
            else
                index_local = nd_index_patch - obj.borderIndex1*ones(1,size(nd_index_patch,2))+1;
                total_size = (obj.borderIndex2 - obj.borderIndex1+1)';
            end
            
            [~,out1] = find(index_local<=0);
            [~,out2]= find(repmat(obj.size,1,size(index_local,2))-index_local<0);
            into_account = setdiff(1:size(index_local,2),unique([out1(:)'  out2(:)']));
            
            if obj.ndimensions==1
                tmp = index_local(1,into_account)';
            else
                str = [];
                for ii=1:ndims(obj.data)
                    str = [str 'index_local(' num2str(ii) ',into_account)'','];
                end
                %[~,out1] = find(index_local - obj.index*ones(1,size(index_local,2))<=0);
                eval(['tmp = sub2ind(total_size,' str(1:end-1) ');' ]);
            end
            
            onedindex = NaN(size(index_local,2),1);
            onedindex(into_account) = tmp;
            
        end
        
        
        function SetOriginToCenter(obj)
            % TODO apply orientation here!
            obj.origin = obj.orientation * (-(obj.size -1).*obj.spacing/2.0 );
        end
        
        function c = GetGeometricalCenter(obj)
            b = obj.GetBounds();
            c = (b(2:2:end)+b(1:2:end))/2;
        end
        
        function out = extractFrame(obj,nframe)
            if numel(obj.size)~=4
                disp('WARNING: input image is not 4D. There might be problems')
            end
            if isa(obj,'PatchType')
                out = PatchType(obj.size(1:end-1),obj.origin(1:end-1),obj.spacing(1:end-1),obj.orientation(1:end-1,1:end-1));
            else
                out = ImageType(obj.size(1:end-1),obj.origin(1:end-1),obj.spacing(1:end-1),obj.orientation(1:end-1,1:end-1));
            end
            out.data = obj.data(:,:,:,nframe);
            
        end
        
        
        
        function h = show(obj,varargin)
            pos = obj.GetPosition(1:prod(obj.size))';
            opacity=-1;
            i=1;
            planeToPlot = []; % three orthogonal planes by default
            crange = [min(obj.data(:)) max(obj.data(:))];
            while (i <= size(varargin,2))
                if  (strcmp( varargin{i} , 'opacity'))
                    opacity = varargin{i+1};
                    i = i+1;
                elseif  (strcmp( varargin{i} , 'colorrange'))
                    crange  = varargin{i+1};
                    i = i+1;
                elseif  (strcmp( varargin{i} , 'plane'))
                    planeToPlot  = varargin{i+1};
                    i = i+1;
                end
                i = i+1;
            end
            
            if numel(obj.size)==2
                X = reshape(pos(:,1),size(obj.data));
                Y = reshape(pos(:,2),size(obj.data));
                Z = zeros(size(obj.data));
                
                data_to_plot=obj.data;
                
                if opacity>=0
                    h= surf(X, Y, Z, data_to_plot, 'FaceColor', 'texturemap','EdgeColor','none','FaceAlpha',opacity);
                else
                    h = surf(X, Y, Z, data_to_plot, 'FaceColor', 'texturemap','EdgeColor','none');
                end
            elseif obj.ndimensions==3
                
                if obj.size(3)==1 % this is a slice of a 3D image
                    %bounds = [min(pos) max(pos)];
                    X = reshape(pos(:,1),size(obj.data));
                    Y = reshape(pos(:,2),size(obj.data));
                    Z = reshape(pos(:,3),size(obj.data));
                    if opacity>=0
                        h = surf(X, Y, Z, obj.data, 'FaceColor', 'texturemap','EdgeColor','none','FaceAlpha',opacity);
                    else
                        h= surf(X, Y, Z, obj.data, 'FaceColor', 'texturemap','EdgeColor','none');
                    end
                else
                    
                    bds = obj.GetBounds();
                    centroid = (bds([1 3 5])+bds([2 4 6]))/2;
                    
                    n(1) = resliceImage(obj,'plane',[1 0 0]',centroid);
                    n(2) = resliceImage(obj,'plane',[0 1 0]',centroid);
                    n(3) = resliceImage(obj,'plane',[0 0 1]',centroid);
                    
                    hold on;
                    if ~numel(planeToPlot)
                        h(1)=n(1).show();
                        h(2)=n(2).show();
                        h(3)=n(3).show();
                    else
                        for ii=planeToPlot
                            h(ii)=n(ii).show();
                        end
                    end
                    hold off;
                    
                    
                end
            else
                disp('show method not implemented for 4D images yet')
            end
            caxis(crange.*[1 1.05]);
        end
        
    end
    
    methods(Access = private)
        function  value = evaluate_at_point_linear(obj,continuous_index)
            % continuous_index is 3 x N
            
            ind = find(obj.size<2);
            
            
            str = '';
            for i = 1:size(continuous_index,1)
                str = [str ', continuous_index('  num2str(i) ',:)'];
            end
            
            
            
            if ~numel(ind)
                % there are no singleton dimensions
                value = eval(['interpn(double(obj.data)' str ',''linear'')']);
            else
                if numel(ind)==1
                    data2 = cat(ind,obj.data,obj.data,obj.data);
                    continuous_index(ind,:)=continuous_index(ind,:)+1;
                    value = eval(['interpn(double(data2)' str ',''linear'')']);
                end
            end
            
            value(isnan(value))=obj.paddingValue;
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

