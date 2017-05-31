function [imout]=normimage(imagein,range,howarg);
% NORMIMAGE - normalize image range
% [imout]=normimage(imagein,range,howarg);
%  range: range to scale image
%  howarg : how to scale, 'scale' (default) or 'wrap'

if nargin < 3,
  howarg = 'scale';
end;

if strcmp(lower(howarg),'scale'),
  % [0 to 1]
  imout = imagein - min(imagein(:));
  imout = imout ./  max(imout(:));

  % range
  if nargin > 1,
    imout = imout .* (max(range)-min(range));
    imout = imout +  min(range);
  end;
else,
  rl = (max(range)-min(range));
  imout = imagein;
  ii = find(imout>max(range));
  while ~isempty(ii),
    imout(ii) = imout(ii)-rl;
    ii = find(imout>max(range));
  end;
  ii = find(imout<min(range));
  while ~isempty(ii),
    imout(ii) = imout(ii)+rl;
    ii = find(imout<min(range));
  end;
end;
