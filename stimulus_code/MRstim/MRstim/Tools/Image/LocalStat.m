function [LocalMean,LocalSd,AveSd]=LocalStat(im,scale)
% LocalStat - compute local rms
%  [LocalMean,LocalSd,AveSd]=LocalStat(im,scale)
%   scale is standard deviation of Gaussian weighting function

% from Steven Dakin

fim=fft2(im);
fim2=fft2(im.^2);
LocalMean=DoGaussian(fim,scale);
LocalMeanSq=DoGaussian(fim2,scale);
for ii=1:length(scale)
        LocalSd(:,:,ii)=sqrt(LocalMeanSq(:,:,ii)-(LocalMean(:,:,ii).^2));
end
AveSd=sqrt(sum(LocalSd(:).^2));
