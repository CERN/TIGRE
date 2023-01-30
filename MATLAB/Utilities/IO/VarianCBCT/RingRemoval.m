function proj = RingRemoval(proj, varargin)
% column-order based filtering
% Method: median or trimmean
% Input:
%       proj: projection matrix
%       varargin: filter type, filter kernel size
% Output:
%       proj: filtered projection matrix
% Date: 2021-09-03
% Author: Yi Du, yi.du@hotmail.com

if(isempty(varargin))
    tag_filter = 'median';
    nkernel = 9;
elseif( length(varargin)==1 )
    if(ischar(varargin{1}))
        if(contains(varargin{1},'median'))
            tag_filter = 'median';
        elseif((contains(varargin{1},'trimmean')))
            tag_filter = 'trimmean';
        else
            error('filter type error');
        end
    else
        error('filter type is missing');
    end
    nkernel = 9;
elseif( length(varargin)== 2)
    if(ischar(varargin{1}))
        if(contains(varargin{1},'median'))
            tag_filter = 'median';
        elseif((contains(varargin{1},'trimmean')))
            tag_filter = 'trimmean';
        else
            error('filter type error');
        end
    else
        error('filter type is missing');
    end
    
    if(isnumeric(varargin{2}))
        nkernel = varargin{2};
    end
    
    if(nkernel ~= round(varargin{2}))
        error('filter kernel size error');
    end        
else
    error('input paramter error');
end


% 3D projection matrix
if(contains(tag_filter,'median'))
    idx = round( (nkernel + 1)/2);
    for ii = 1:size(proj,3)
        proj(:,:,ii) = ordfilt2(proj(:,:,ii), idx, ones(1,nkernel));
    end
elseif(contains(tag_filter,'trimmean'))
    nb = round( (nkernel+1)/2) -1;
    left = nb+1;
    right = size(proj, 2) - left;
    for ii = left: right
        tmp_blk = proj(:, ii-nb : ii+nb, :);
        proj(:,ii,:) = trimmean(tmp_blk, 40, 2, 'round');
    end
end
    
end

