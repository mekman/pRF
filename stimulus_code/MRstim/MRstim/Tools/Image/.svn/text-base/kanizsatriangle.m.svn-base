function k=kanizsatriangle(sz,w)
% kanizsatriangle - make Kanizsa triangle
%
% k=kanizsatriangle(sz,w)
%
% input:
% sz - size of image
% w - size of weighted cosine window (0 = none)
%
% 2010/01 SOD: wrote it.

if ~exist('sz','var') || isempty(sz)
    sz = 512;
end
if ~exist('w','var') || isempty(w)
    w = false;
end


% make image
k = zeros(sz); % sz = diameter

% make circles
c = makecircle(round(sz/6));
szc= size(c,1)/2;

% Make matrix coordinate system with 0 in center:
[xx,yy]=meshgrid(-sz/2:1:sz/2);
[th r]=cart2pol(xx,yy);

% find angles to place circles
angles = linspace(0,2*pi,4);
angles = angles(1:3)-pi/2;

% compute xy-coordinates
[x,y] = pol2cart(angles,sz/4*[1 1 1]);
x = x + sz./2;
y = y + sz./2;


% place circles
for n=1:numel(angles)
    cx = round([x(n)-szc, x(n)+szc]);
    cy = round([y(n)-szc, y(n)+szc]);
    k(cy(1)+1:cy(2),cx(1)+1:cx(2))=c;
end

% mask triangle
mask = poly2mask(x,y,sz,sz);
k(mask) = 0;
imagesc(mask)

% create cosine window
if w>0
    % center around 0
    k = -k-0.5;
    % make window
    win = makecircle(sz-2*w,sz,w);
    % get weighted stats
    s = wstat(k(:),win(:));
    % center on mean window
    k = k - s.mean;
    % scale to max contrast
    k = k./max(abs(k(win>0)));
    % mask
    k = k.*win;
    % scale to 0-1
    k = k./2+0.5;
    % convert to uint8
    k = uint8(k.*255);
end







