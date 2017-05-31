function [imout] = simulateCMF2pRF;
% fromimagetov1 - warp image to approx v1 size
% Uses log-transform, see articles by Schwartz et al.
%

% default params
p = sub_default_params;

% make visual image
im = sub_visual_field_image(params)

% transform visual image to cortical representaton
[im inv1] = fromimagetov1(im);

% make cortical RFs that are a constant size on the cortex (Gaussian or
% circular)
rfs = % 2d cRFimage x number pixels in V1

% transform cortical RFs to visual field
rfs = fromv1toimage(rfs);

% compute visual field sizes
rfs = computeRFsize(rfs)


% provide output
imout.ecc = im.ecc;
imout.pol = im.pol;
imout.pRF = rfs(inv1);
imout.mm  = pixelcoordinates;
imout.roi 
return


function sub_default_params
% some parameters
imsize.pix = 1024;
imsize.deg = 20;
lp.a = 1.5; % foveal starting angle
lp.b = 30;  % eccentric tapering
lp.k = 0.01;
lp.s = 4; % eccentricity scaling
return

function sub_visual_field_image(p)
% make visual field image

scale2deg = (imsize.deg.*2)./imsize.pix;
x = linspace(-imsize.pix./2.*scale2deg,imsize.pix./2.*scale2deg,imsize.pix);
[X Y] = meshgrid(fliplr(x),x);
return



% convert VF image to polar coordinates
Z = X+Y.*sqrt(-1);
th = angle(Z);
r  = abs(Z);
figure(2);
subplot(2,3,1);imagesc(mod(th,2*pi)-pi);axis image off;colorbar
subplot(2,3,2);imagesc(r,[0 imsize.deg]);axis image off;colorbar

% log conversion from VF to V1
a=1.5 % foveal starting angle
b=30  % eccentric tapering
k=0.01
s=4 % eccentricity scaling
W = k*log((Z+a)./(Z+b));
minZ = min(r(:));
minW = k*(log((minZ+a)./(minZ+b)));
maxW = k*(log((imInDeg+a)./(imInDeg+b)));
W2 = (W-minW)./(maxW-minW).*imInDeg;

% original log scaling
%a=1.5
%W = (k*log((Z+a)));
%minW = (k*(log((0+a))))
%maxW = (k*(log((20+a))))
%W = (W-minW)./(maxW-minW).*20;


% get new polar coordinates, a bit of black magic here
th2 = angle(W2-Z);
r2 = abs(W2-Z);
[newr, Y2] = pol2cart(th2,r2);
newth = Y2./max(Y2(:))*pi; 

% scale ecc
mask1 = abs(mod(newth,2*pi)-pi)<=pi./2;
mask1 = mask1(:,round(sz2./2):end);
rmask = newr(:,round(sz2./2):end);
minW = min(rmask(mask1));
newr = newr-minW;
rmask = rmask-minW;
maxW = max(rmask(mask1));
pos  = newr>0;
newr(pos) = ((newr(pos)./maxW).^s).*(imInDeg+1);

figure(2);
subplot(2,3,4);imagesc(mod(newth,2*pi)-pi);axis image off;colorbar
subplot(2,3,5);imagesc(newr,[0 imInDeg]);axis image off;colorbar

% now make pixel coordinates
x = 1:sz2;
y = 1:sz1;
[X Y] = meshgrid(x,y);%figure(3);imagesc(X);colorbar
maskold = abs(mod(th,2*pi)-pi)<=pi./2 & r>=0 & r<=imInDeg;
masknew = abs(mod(newth,2*pi)-pi)<=pi./2 & newr>=0 & newr<=imInDeg;
figure(2);
subplot(2,3,3);imagesc(maskold);axis image off;colorbar
subplot(2,3,6);imagesc(masknew);axis image off;colorbar

% reshape only values within mask
th2 = newth(masknew);
r2  = newr(masknew);

% nearest neighbor interpolation
[x,y] = pol2cart(th2,r2);
x = round((-x./imInDeg./2+.5)*(sz2-1)+1);
y = round((-y./imInDeg./2+.5)*(sz1-1)+1);
imout = zeros(prod(size(img)),1)+bgindex;
imout(masknew) = img(sub2ind([sz1 sz2],y,x));
imout = reshape(imout, sz1, sz2);   
orgim = double(img).*maskold;
orgim(~maskold) = bgindex;

% plot results
figure(1);
subplot(1,2,1);imagesc(orgim);colormap gray;axis image off;
subplot(1,2,2);imagesc(imout);colormap gray;axis image off;
    
return
