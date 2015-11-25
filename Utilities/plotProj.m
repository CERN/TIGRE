function plotProj(proj,alpha,varargin)
% PLOTPROJ plots the projection data
% 
%   PLOTPROJ(PROJ,ALPHA) plots projection data PROJ over radians ALPHA,
%                        where length(ALPHA) == size(PROJ,3)
% 
%   PLOTPROJ(PROJ,ALPHA,OPTS,VAL,...) Allows user specified OPTS options
%                                     with corresponding VAL values.
%           Posible options in OPTS are:
% 
%   'Colormap': Sets the colormap. Possible values for VAL are the names of
%               the stadard MATLAB colormaps, the names in the perceptually 
%               uniform colormaps tool or a custom colormap, being this last 
%               one a 3xN matrix. Default is GRAY
%   'Clims':    a 2x1 matrix setting the upper and lower limits of the
%               colors. The default is the limits of the data.
%   'Step':     Sets the step size between slice and slice. Step is 1 by
%               default.
%   'Savegif':  With an string in VAL, saves the image as .gif with
%               VAL as filename
%
%
%% Parse inputs
opts=     {'Step','Colormap','Clims','Savegif'};
defaults= [    1,    1     ,   1 ,      1];

% Check inputs
nVarargs = length(varargin);
if mod(nVarargs,2)
    error('CBCT:plotProj:InvalidInput','Invalid number of inputs')
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
        val=varargin{ninput*2};
        ninput=ninput+1;
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
                steps=varargin{ii+1};
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
                climits=[min(proj(:)) max(proj(:))];
            else
                if min(size(val))==1 && max(size(val))==2
                    climits=val;
                else
                    error('CBCT:plotImgs:InvalidInput','Invalid size of Clims')
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
        otherwise
          error('CBCT:plotImgs:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in plotImg()']);
    end
end

%% Do ploting
fh=figure();

for ii=1:steps:size(proj,3)
    image=squeeze(proj(:,:,ii));
    imagesc((image));
    
    axis image;  
    axis equal; 
    
    colormap(cmap); 
    colorbar; 
    caxis([climits(1),climits(2)]);
    
    xlabel('-> U');
    ylabel('-> V');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'YDir','normal');
    
    title(['Degree : ',num2str(alpha(ii)*180/pi)]);
    pause(0.01);
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
