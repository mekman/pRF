function [imout, orgim, orgimtot] = fromimagetov1(img,bgindex);
% fromimagetov1 - warp image to approx v1 size
% Uses log-transform, see articles by Schwartz et al.
% 
% out = fromimagetov1(img,bgindex);
%
% 2007/01 SOD: wrote it.

if nargin < 1, % for debugging
    img = imread('eight.tif')+1;
end;
if nargin < 2,
    bgindex =0;
end;

origimg = img;
% upsample and size
img = imresize(origimg,2,'bicubic');
[sz1, sz2] = size(img);

% make cartesian coordinates
% assume the image covers first 120 degrees radius of the visual
% field
imInDeg = 20;
scale2deg = (imInDeg.*2)./min(sz1,sz2);
x = linspace(-sz2./2.*scale2deg,sz2./2.*scale2deg,sz2);
y = linspace(-sz1./2.*scale2deg,sz1./2.*scale2deg,sz1);
[X Y] = meshgrid(fliplr(x),y);%figure(3);imagesc(X);colorbar

% convert to polar coordinates
Z = X+Y.*sqrt(-1);
th = angle(Z);
r  = abs(Z);
figure(2);
subplot(2,3,1);imagesc(mod(th,2*pi)-pi);axis image off;colorbar
subplot(2,3,2);imagesc(r,[0 imInDeg]);axis image off;colorbar

if nargout>2,
    ii = r>imInDeg;
    orgimtot = img;
    orgimtot(ii) = bgindex;
end;

% log conversion
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
