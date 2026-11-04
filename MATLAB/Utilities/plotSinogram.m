function plotSinogram(proj,alpha,varargin)
% PLOTSINOGRAM plots the sinogram of the projection data
%
% PLOTSINOGRAM(PROJ,ALPHA) plots the sinogram of projection data PROJ
% over radians ALPHA, where length(ALPHA) == size(PROJ,3). This does the
% same as PLOTPROJ, but with the input transposed and sliced along a
% different axis: instead of showing, for each angle, the U-V detector
% image, PLOTSINOGRAM shows, for each detector row (V), the U-angle
% sinogram (i.e. the detector's horizontal readout stacked over all
% acquisition angles).
%
% PLOTSINOGRAM(PROJ,ALPHA,OPTS,VAL,...) Allows user specified OPTS options
% with corresponding VAL values.
% Posible options in OPTS are:
%
% 'Colormap': Sets the colormap. Possible values for VAL are the names of
% the standard MATLAB colormaps, the names in the perceptually
% uniform colormaps tool or a custom colormap, being this last
% one a 3xN matrix. Default is GRAY
% 'Clims': a 2x1 matrix setting the upper and lower limits of the
% colors. The default is the limits of the data.
% 'Step': Sets the step size between slice and slice. Step is 1 by
% default.
% 'Savegif': With an string in VAL, saves the image as .gif with
% VAL as filename
% 'Slice' : Plot a single sinogram (single detector row), given by
% index VAL. Should not be used with 'Step', overwrites it.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
%
% Copyright (c) 2015, University of Bath and
% CERN-European Organization for Nuclear Research
% All rights reserved.
%
% License: Open Source under BSD.
% See the full license at
% https://github.com/CERN/TIGRE/blob/master/LICENSE
%
% Contact: tigre.toolbox@gmail.com
% Codes: https://github.com/CERN/TIGRE/
% Coded by: Ander Biguri
%--------------------------------------------------------------------------
%% Parse inputs
opts= {'step','colormap','clims','savegif','slice'};
defaults= [ 1, 1 , 1 , 1 , 1];

% Check inputs
nVarargs = length(varargin);
if mod(nVarargs,2)
    error('CBCT:plotSinogram:InvalidInput','Invalid number of inputs')
end
% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(varargin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    else
        error('CBCT:plotSinogram:InvalidInput',['Optional parameter "' varargin{ii} '" does not exist' ]);
    end
end

for ii=1:length(opts)
    opt=opts{ii};
    default=defaults(ii);
    % if one option isnot default, then extranc value from input
    if default==0
        ind=[];jj=1;
        while isempty(ind)
            ind=find(isequal(opt,lower(varargin{jj})));
            jj=jj+1;
        end
        if isempty(ind)
            error('CBCT:plotSinogram:InvalidInput',['Optional parameter "' varargin{jj} '" does not exist' ]);
        end
        val=varargin{jj};
    end
    switch opt
        % % % % % % %Step
        case 'step'
            if default
                steps=1;
            else
                if ~isnumeric(val)
                    error('CBCT:plotSinogram:InvalidInput','Invalid step')
                end
                steps=val;
            end
            % % % % % % % Colormap choice
        case 'colormap'
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
                        error('CBCT:plotSinogram:InvalidInput','Invalid size of colormap')
                    end
                    cmap=val;
                end
            end
            % % % % % % % Limits of the colors
        case 'clims'
            if default
                climits=[min(proj(:)) max(proj(:))];
            else
                if min(size(val))==1 && max(size(val))==2
                    climits=val;
                else
                    error('CBCT:plotSinogram:InvalidInput','Invalid size of Clims')
                end
            end
            % % % % % % % % do you want to save result as gif?
        case 'savegif'
            if default
                savegif=0;
            else
                savegif=1;
                if ~ischar(val)
                    error('CBCT:plotSinogram:InvalidInput','filename is not character')
                end
                filename=val;
            end
        case 'slice'
            if default
                slice=0;
            else
                slice=val;
            end
        otherwise
            error('CBCT:plotSinogram:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in plotSinogram()']);
    end
end

%% Do ploting
fh=figure();

if slice
    list=slice;
else
    list=1:steps:size(proj,1);
end

for ii=list
    % Transpose relative to plotProj: fix a detector row (dim 1) and
    % show the U (dim 2) vs angle (dim 3) sinogram for that row.
    image=squeeze(proj(ii,:,:));
    imagesc((image));

    axis image;
    axis equal;

    colormap(cmap);
    colorbar;
    caxis([climits(1),climits(2)]);

    xlabel('-> Angle');
    ylabel('-> U');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'YDir','normal');

    title(['Detector row : ',num2str(ii)]);
    drawnow;
    pause(0.01);
    if savegif

        frame = getframe(fh);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if ii == list(1)
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
        end
    end
end
