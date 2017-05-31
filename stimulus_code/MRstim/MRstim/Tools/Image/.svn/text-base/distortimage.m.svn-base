function imout = distortimage(im,dispx,dispy);
% DISTORTIMAGE - distort image according to displacement fields (x,y)
% imout = distortimage(im,dispx,dispy);

% adjusted from Ethan Meyers 2002  emeyers@mit.edu 

if nargin < 3,
  help(mfilename);
  return;
end;

[m,n]=size(im);

[coordsX coordsY] = meshgrid(1:n, 1:m);

newCoordsX = coordsX + dispx;
newCoordsY = coordsY + dispy;

% adjust any pixels that might be out of range of the image
newCoordsX(find(newCoordsX < 1)) = 1; 
newCoordsX(find(newCoordsX > n)) = n; 
newCoordsY(find(newCoordsY < 1)) = 1; 
newCoordsY(find(newCoordsY > m)) = m; 

linearX = reshape(newCoordsX, 1, m*n); 
linearY = reshape(newCoordsY, 1, m*n); 

indecies = sub2ind([m,n], linearY, linearX);  

im = reshape(im, 1, n*m); % 1d

imout = im(indecies);
imout = reshape(imout, m, n);   % reshape back into an image

