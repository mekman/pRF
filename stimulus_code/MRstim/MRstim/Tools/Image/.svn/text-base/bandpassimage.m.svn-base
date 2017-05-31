function [imout]=bandpassimage(imagein,freqrange,scale);
% BANDPASSIMAGE - bandpass filter an image

if nargin<3, % scale 0-1 input to imagerange
  freqrange = freqrange .* size(imagein,1)/2
end;

padd=size(imagein,1)*3;

% circular coords
coords=ceil(-padd/2:padd/2-1);
[xx,yy]=meshgrid(coords);
dist_from_center = round(sqrt((xx.*xx)+(yy.*yy)));
 
meanf = mean(imagein(:));
fftimagein = fftshift( fft2(imagein-meanf,padd,padd) );
imout = zeros(size(imagein,1),size(imagein,1),size(freqrange,1));

for n=1:size(freqrange,1),
  % make filters
  filt = zeros(padd);
  filt(dist_from_center>=freqrange(n,1) & dist_from_center<=freqrange(n,2))=1;

  % filter
  imouttmp = real(ifft2(fftshift( fftimagein .*filt ),padd,padd))+meanf;
  imout(:,:,n) = imouttmp(1:size(imagein,1),1:size(imagein,2));
end
