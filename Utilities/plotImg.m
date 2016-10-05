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
%   'Slice'     Plots a single slice of the data. Overwrittes Step.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
% 
% Copyright (c) 2015, University of Bath and 
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License:            Open Source under BSD. 
%                     See the full license at
%                     https://github.com/CERN/TIGRE/license.txt
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Ander Biguri
%--------------------------------------------------------------------------
%% Parse inputs

opts=     {'Step','Dim','Savegif','Colormap','Clims','Slice'};
defaults= [   1  ,  1  ,    1    ,    1     ,   1 ,1];

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

ninput=1;
for ii=1:length(opts)
    opt=opts{ii};
    default=defaults(ii);
    % if one option isnot default, then extranc value from input
    if default==0
        ind=double.empty(0,1);jj=1;
        while isempty(ind)
            ind=find(isequal(opt,varargin{jj}));
            jj=jj+1;
        end
        val=varargin{jj};
    end
    
    switch opt
% % % % % %         %Step 
        case 'Step'
            if default
                steps=1;
            else
                if ~isnumeric(val)
                    error('CBCT:plotImgs:InvalidInput','Invalid step')
                end
                steps=val;
            end
% % % % % %         % Plot single slice
         case 'Slice'
            if default
                slice=0;
            else
                slice=val;
            end

% % % % % %         % iterate trhoug what dim?
        case 'Dim'
            if default
                crossect=1;
            else
                if length(val)>1  
                    error('CBCT:plotImgs:InvalidInput','Invalid Dim')
                end
                
                if val==3 || lower(val)=='z'
                    crossect=3;
                end
                if val==2 || lower(val)=='y'
                    crossect=2;
                end
                if val==1 || lower(val)=='x'
                    crossect=1;
                end
            end
            
% % % % % % %         % do you want to save result as gif?
        case 'Savegif'
            if default
                savegif=0;
            else
               savegif=1;
               if ~ischar(val)
                   error('CBCT:plotImgs:InvalidInput','filename is not character')
               end
               filename=val;
            end
% % % % % %         % Colormap choice
        case 'Colormap'
            if default
                cmap='gray';
            else
                
                if ~isnumeric(val)  
                    % check if it is from perceptually uniform colormaps.
                    if ismember(val,{'magma','viridis','plasma','inferno'})
                        cmap=eval([val,'()']);
                    else
                        cmap=val;
                    end                   
                else
                    % if it is a custom colormap
                    if size(val,2)~=3
                        error('CBCT:plotImgs:InvalidInput','Invalid size of colormap')
                    end
                    cmap=val;
                end
            end
% % % % % %         % Limits of the colors
        case 'Clims'
            if default
                if any(size(img)>512) || ~areTheseToolboxesInstalled({'MATLAB','Statistics Toolbox'})  || areTheseToolboxesInstalled({'MATLAB','Statistics and Machine Learning Toolbox'})
                   climits=[0, 0.95*max(img(:))];
                else
                   climits=prctile(double(img(:)),[1 99]);
                end
            else
                if min(size(val))==1 && max(size(val))==2
                    climits=val;
                else
                    error('CBCT:plotImgs:InvalidInput','Invalid size of Clims')
                end
            end
            
        otherwise
          error('CBCT:plotImgs:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in plotImg()']);
    end
end



%% Do the plotting!

fh=figure();

if slice
    list=slice;
else 
    list=1:steps:size(img,crossect);
end
for ii=list
    
    if crossect==1
      imagesc(squeeze(img(ii,:,:)).'); 
    else if crossect==2
             imagesc(squeeze(img(:,ii,:)).');
        else
             imagesc(squeeze(img(:,:,ii)));
        end
    end
    axis image; 
    axis equal; 
    
    colormap(cmap); 
    colorbar; 
    caxis([climits(1),climits(2)]);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'YDir','normal');
    
    
    if crossect==3 
        xlabel('->Y');
        ylabel('X<-');
        title(['Bottom to top->Z : ',num2str(ii)]);
    end
     if crossect==2 
        xlabel('->X');
        ylabel('->Z');
        title(['Left to Rigth direction ->Y : ',num2str(ii)]);
    end
    if crossect==1 
        xlabel('->Y');
        ylabel('->Z');
        title(['Detector to Source direction ->X : ',num2str(ii)]);
    end
    drawnow
    
    if savegif
        
      frame = getframe(fh);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if ii == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
      end
    end
end

