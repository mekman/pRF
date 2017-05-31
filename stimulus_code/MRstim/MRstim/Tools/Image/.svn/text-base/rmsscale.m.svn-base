function [scalefac,scalematrix]=rmsscale(matrix,rms,range);
% RMSSCALE - scale image to rms-contrast
%  [scalefac,scalematrix]=rmsscale(matrix,rms,range);
%  default range [0 255]

if nargin < 2,
  help(mfilename);
  return;
end;
if nargin < 3 | isempty(range),
  range = [0 255];
end;

hrange=(range(2)-range(1))./2;
disc=imagestat(matrix);
out=matrix;
scalefac=1;
scalefac2=0;

while round((disc.contrast_rms-rms).*1000)~=0,
  scalefac2=disc.contrast_rms.\rms;
  scalefac=scalefac.*scalefac2;
  
  meanmatrix=mean(matrix(:));
  out = (matrix-meanmatrix).*scalefac+meanmatrix;
  out(find(out>range(2)))=range(2);
  out(find(out<range(1)))=range(1);
  disc=imagestat(out);
  disc.contrast_rms;
end;
scalematrix=out;
 
  
