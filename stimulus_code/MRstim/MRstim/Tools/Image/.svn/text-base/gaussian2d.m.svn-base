function [gauss]=gaussian2d(x,y,sigmax,sigmay,theta,derivative);
% GAUSSIAN2D - 2D Gaussian and/or derivatives
% gauss=gaussian2d(x,y,sigmax,sigmay,theta,derivative);

if nargin < 3,
    help(mfilename);
    return;
end;
if nargin < 4 | isempty(sigmay),
    sigmay = sigmax;
end;
if nargin < 5 | isempty(theta),
  theta = 0;
end;
if nargin < 6 | isempty(derivative),
    derivative = [];
end;

% make grid
xx=1:x;xx=xx-mean(xx);
yy=1:y;yy=yy-mean(yy);
[xtmp,ytmp]=meshgrid(xx,yy);
gx = xtmp.*cos(theta) + ytmp.*sin(theta);
gy = ytmp.*cos(theta) - xtmp.*sin(theta);

% make gaussian
Gx = exp(-0.5*gx.^2/sigmax^2);
Gy = exp(-0.5*gy.^2/sigmay^2);
gauss = Gx.*Gy;
gauss = gauss./sum(gauss(:));


% make derivatives
if ~isempty(derivative),
    gauss.gauss = gauss;
    gauss.dx    = -gx.*gauss.gauss./sigmax.^2;
    gauss.dy    = -gy.*gauss.gauss./sigmay.^2;
end; 


