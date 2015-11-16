%% Colormaps that are pretty awesome: DEMO2
%  In this demo we will show how to supereasily use the new colormaps
%
%
% this is an exampe of 4 colormaps that are being considered as the default
% colormap in python's matplotlib lybrary.
%
% All of them look quite good and they dont have any official name, so at
% the moment they are A,B,C,D.
%
% colormaps from https://github.com/bids/colormap
%
% Ander Biguri
%% Clear workspace and get screen data
clear;
clc
close all;


%% Generate sample data
X=peaks(200);

%% Load Colomaps

jet=colormap('jet');
parula=parula();
magma=magma();
inferno=inferno();
plasma=plasma();
viridis=viridis();



%% Chose colormap
% Use 1 only, else it will  just use the last
% CTRL+R -> comment line
% CTRL+T -> Uncomment line

% colormap(jet);
% colormap(parula);
% colormap(magma);
% colormap(inferno);
% colormap(plasma);
colormap(viridis);

%% Plot

for ii=1:0.3:20
    surf(cos(2*pi*ii/20)*X,'linestyle','none');
    axis off
    axis([0 200 0 200 -10 10])
    set(gca, 'CLim', [-8, 8]);
    drawnow
end
