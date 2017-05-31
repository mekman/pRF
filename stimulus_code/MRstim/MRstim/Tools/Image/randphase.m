function [out,pow] = randphase(matrix,rand_orient,padd,mask,nphaseshifts,doPlot);
% RANDPHASE - phase (and orientation) scrambling
% [out] = randphase(matrix,rand_orient,padd,mask,nphaseshifts,doPlot);

% ---------------------------------
% Copyright (c) 2003 Serge Dumoulin
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License version 2
% as published by the Free Software Foundation.
% ---------------------------------

if nargin < 1,
   help(mfilename);
   return;
end;
if nargin < 2 | rand_orient == 0,
   rand_orient = [];
end;
if nargin < 3 | padd == 0,
   padd = size(matrix,1);
end;
if nargin < 4,
   mask = [];
end;
if nargin < 5,
   nphaseshifts = 0;
end;
if nargin < 6 | doPlot == 0,
   doPlot = [];
end;
rand('state',sum(100*clock));

[len,wid]=size(matrix);
if len ~= wid,
   disp(sprintf('Error --- not a square matrix'));
end;

% get fourierspectrum
mean_matrix=mean(matrix(:));
fft_matrix=matrix- mean_matrix;
if nphaseshifts==0,
  if ~isempty(mask),
    fft_matrix = fft_matrix .* mask;
  end;
  fft_matrix=fftshift(fft2(fft_matrix,padd,padd));
  pow=abs(fft_matrix);
  old_fft_matrix=fft_matrix;
else,
  ps = zeros(padd,padd,(nphaseshifts+1));
  range = max(fft_matrix(:))-min(fft_matrix(:));
  shift = range/(nphaseshifts+1);
  for n=1:(nphaseshifts+1),
    t0 = normimage(fft_matrix + (n-1)*shift,...
                   [min(fft_matrix(:)) max(fft_matrix(:))],...
                   'wrap');
    if ~isempty(mask),
      t0 = t0 .* mask;
    end;

    %figure(1);clf;imagesc(t0.*mask);axis('image','off');
    %title(sprintf('%d',n));colorbar;drawnow;
    
    tmp = fftshift(fft2(t0,padd,padd));
    ps(:,:,n)  = abs(tmp); 
  end;
  pow=median(ps,3);
  old_fft_matrix=fft_matrix;
end;

if ~isempty(rand_orient),
   new_pow=pow;
   % circular coordinate system
   coords=ceil(-padd/2:padd/2-1);
   [xx,yy]=meshgrid(coords);
   dist_from_center = sqrt((xx.*xx)+(yy.*yy));
   % distribute energy in circles
   for n=0:round(max(dist_from_center(:))),
      target_circle=find(dist_from_center >= n-.5 & ...
                        dist_from_center <  n+.5 );
      if length(target_circle) == 1,
        new_pow(target_circle)=mean(pow(target_circle));
      else,
        new_pow(target_circle)=mean(pow(target_circle)).*...
            repmat(rand(length(target_circle)/2,1)+.5,2,1);
      end;
   end;
   pow=new_pow;
end;

rand_ang=rand(padd).*(2*pi); 
newmatrix=pow.*(cos(rand_ang)+sin(rand_ang).*sqrt(-1));

out=real(ifft2(ifftshift(newmatrix),padd,padd))+mean_matrix;
out=out(1:len,1:wid);

if ~isempty(doPlot),
figure;
set(gcf,'Position',[100 100 500 800],'Color',[1 1 1]);
subplot(5,2,1);imagesc(matrix(1:len,1:wid));
title('org - image');axis('square');colormap('gray');colorbar;
subplot(5,2,2);imagesc(out(1:len,1:wid));
title('new - image');axis('square');colormap('gray');colorbar;
subplot(5,2,3);imagesc(abs(old_fft_matrix(1:len,1:wid)));
title('org - abs');axis('square');colormap('gray');colorbar;
subplot(5,2,4);imagesc(abs(newmatrix(1:len,1:wid)));
title('new - abs');axis('square');colormap('gray');colorbar;
subplot(5,2,5);imagesc(angle(old_fft_matrix(1:len,1:wid)));
title('org - angle');axis('square');colormap('gray');colorbar;
subplot(5,2,6);imagesc(angle(newmatrix(1:len,1:wid)));
title('new - angle');axis('square');colormap('gray');colorbar;
subplot(5,2,7);imagesc(real(old_fft_matrix(1:len,1:wid)));
title('org - real');axis('square');colormap('gray');colorbar;
subplot(5,2,8);imagesc(real(newmatrix(1:len,1:wid)));
title('new - real');axis('square');colormap('gray');colorbar;
subplot(5,2,9);imagesc(imag(old_fft_matrix(1:len,1:wid)));
title('org - imag');axis('square');colormap('gray');colorbar;
subplot(5,2,10);imagesc(imag(newmatrix(1:len,1:wid)));
title('org - imag');axis('square');colormap('gray');colorbar;
end;

