function [im] = jointhist(x,y,z);
% JOINTHIST - Joint histogram of two images
% [im] = jointhist(y,x,z);

if nargin <2
  help(mfilename);
  return;
end
if nargin < 3
    z = 128;
end;

y = round(normimage(y,[1 z],'scale'));
x = round(normimage(x,[1 z],'scale'));

[rows, cols]=size(y);
im = zeros(z);


% slow loop
fprintf('[%s]Processing:',mfilename);drawnow;
for n1=1:rows,
  for n2=1:cols,
    im(x(n1,n2),y(n1,n2))= im(x(n1,n2),y(n1,n2))+1;
  end;
end;
fprintf('Done.\n');drawnow;
  
if ~nargout,
  imagesc(im);colormap('gray');axis('image');
end;
