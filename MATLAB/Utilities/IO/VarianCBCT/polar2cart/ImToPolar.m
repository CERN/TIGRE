function imP = ImToPolar (imR, rMin, rMax, M, N)
% IMTOPOLAR converts rectangular image to polar form. The output image is 
% an MxN image with M points along the r axis and N points along the theta
% axis. The origin of the image is assumed to be at the center of the given
% image. The image is assumed to be grayscale.
% Bilinear interpolation is used to interpolate between points not exactly
% in the image.
%
% rMin and rMax should be between 0 and 1 and rMin < rMax. r = 0 is the
% center of the image and r = 1 is half the width or height of the image.
%
% V0.1 7 Dec 2007 (Created), Prakash Manandhar pmanandhar@umassd.edu

[Mr Nr] = size(imR); % size of rectangular image
Om = (Mr+1)/2; % co-ordinates of the center of the image
On = (Nr+1)/2;
sx = (Mr-1)/2; % scale factors
sy = (Nr-1)/2;

imP  = zeros(M,  N);

delR = (rMax - rMin)/(M-1);
delT = 2*pi/N;

% loop in radius and 
for ri = 1:M
for ti = 1:N
    r = rMin + (ri - 1)*delR;
    t = (ti - 1)*delT;
    x = r*cos(t);
    y = r*sin(t);
    xR = x*sx + Om;  
    yR = y*sy + On; 
    imP (ri, ti) = interpolate (imR, xR, yR);
end
end

function v = interpolate (imR, xR, yR)
    xf = floor(xR);
    xc = ceil(xR);
    yf = floor(yR);
    yc = ceil(yR);
    if xf == xc & yc == yf
        v = imR (xc, yc);
    elseif xf == xc
        v = imR (xf, yf) + (yR - yf)*(imR (xf, yc) - imR (xf, yf));
    elseif yf == yc
        v = imR (xf, yf) + (xR - xf)*(imR (xc, yf) - imR (xf, yf));
    else
       A = [ xf yf xf*yf 1
             xf yc xf*yc 1
             xc yf xc*yf 1
             xc yc xc*yc 1 ];
       r = [ imR(xf, yf)
             imR(xf, yc)
             imR(xc, yf)
             imR(xc, yc) ];
       a = A\double(r);
       w = [xR yR xR*yR 1];
       v = w*a;
    end
