function plotImg(img,varargin)
% PLOTIMG plots a 3D image in slices
%   PLOTIMG(IMG) plots 3D image IMG looping thourgh X axis (first
%   dimension)
%   PLOTIMG(IMG,OPT,VAL,...) uses options and values for plotting. The
%   possible options in OPT are:
%
%   'Step':     Sets the step size between slice and slice. Step is 1 by
%               default.
%   
%   'Dim':      Sets the dimensions in which the function iterates trhough.
%               Default is 1 or 'X', possibilities are 2,3 or 'Y','Z' 
%               respectively. 
%   
%   'Savegif':  With an string in VAL, saves the image as .gif with
%               VAL as filename
%   'Colormap': Sets the colormap. Possible values for VAL are the names of
%               the stadard MATLAB colormaps, the names in the perceptually 
%               uniform colormaps tool or a custom colormap, being this last 
%               one a 3xN matrix. Default is GRAY
%   'Clims':    a 2x1 matrix setting the upper and lower limits of the
%               colors. The default computes the lower and upper percentile
%               of data, in 1% and 99% and sets the limits to that.

%% Parse inputs

opts=     {'Step','Dim','Savegif','Colormap','Clims'};
defaults= [   1  ,  1  ,    1    ,    1     ,   1 ];

% Check inputs
nVarargs = length(varargin);
if mod(nVarargs,2)
    error('CBCT:plotImgs:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,varargin{ii}));
    if ~isempty(ind)
       defaults(ind)=0; 
    end
end

for ii=1:length(opts)
    opt=otps{ii};
    default=defaults(ii);
    switch opt
        case 'Step'
            if default
            else
            end
        case 'Dim'
            if default
            else
            end
        case 'Savegif'
            if default
            else
            end
            % Colormap
        case 'Colormap'
            if default
                cmap='gray';
            else
                if ismember(varargin{ii+1},{'magma','viridis','plasma','inferno'})
                    cmap=eval('')
                end
            end
            % Limits of the colors
        case 'Clims'
            if default
                climits=prctile(img(:),[1 99]);
            else
                if min(size(varargin{ii+1}))==1 && max(size(varargin{ii+1}))==2
                    climits=varargin{ii+1};
                else
                    error('CBCT:plotImgs:InvalidInput','Invalid size of Clims')
                end
            end
            
        otherwise
           
    end
end



%% Do the plotting!

fh=figure();
for ii=size(img,1):-1*steps:1
    imagesc((squeeze(img(ii,:,:)))'); 
    
    
    axis image; 
    axis equal; 
    
    colormap(cmap); 
    colorbar; 
    caxis([climits(1),climits(2)]);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'YDir','normal');
    
    
    if cross==3 
        xlabel('->Y');
        ylabel('X<-');
        title(['Top to bottom ->Z : ',num2str(ii)]);
    end
    if cross==1 
        xlabel('->Y');
        ylabel('->Z');
        title(['Source to Detector direction ->X : ',num2str(ii)]);
    end
    drawnow
    
    if savegif
        
      frame = getframe(fh);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if ii == 256;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
      end
    end
end

