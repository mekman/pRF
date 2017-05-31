function circularaverage(matrix,corners,doPlot);
% CIRCULARAVERAGE - average in circles


if nargin < 1,
   help(mfilename);
   return;
end;
if nargin < 2 | isempty(corners),
   corners='n';
end;
if nargin < 3 | isempty(doPlot),
   doPlot='n';
end;

% check
[len,wid]=size(matrix);
if len ~= wid,
   disp(sprintf('Error --- not a square matrix'));
end;

% get rel freq 
if isempty(dist_from_center),
  coords=ceil(-len/2:wid/2-1);
  [xx,yy]=meshgrid(coords);
  dist_from_center = round(sqrt((xx.*xx)+(yy.*yy)));
  tmpvar = dist_from_center;
end;

% corners?
if corners == 'y',
   endspec=round(max(dist_from_center(:)));
else,
   endspec=round(len/2);
end;   

