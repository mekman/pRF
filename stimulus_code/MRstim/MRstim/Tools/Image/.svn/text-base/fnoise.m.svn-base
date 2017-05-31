function im=fnoise(Size,f)
% fnoise - make 1/f noise
%
%
% 2009/10 SOD: wrote it.

if ~exist('Size','var') || isempty(Size)
    Size = 512;
end
if ~exist('f','var') || isempty(f)
    f = 2;
end

% make circular coords
coords=-Size/2:Size/2-1;
[xx,yy]=meshgrid(coords);
fractal = sqrt((xx.*xx)+(yy.*yy));
clear xx yy;

fractal(fractal~=0)=1./(fractal(fractal~=0).^f);
fractal = fractal./sum(fractal(:));

% random
im = rand(Size);

% adjust power
im = real(ifft2(fftshift(fftshift(fft2(im)).*fractal)));

if ~nargout
    imagesc(im);colorbar(gray);axis image off;
end

return