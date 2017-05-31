function [powerspectrum,tmpvar] = powerspec(matrix,doPlot,corners,...
    tmpvar,dofft);
% POWERSPEC - gets powerspectrum of (2D) image
% [powerspectrum,tmpvar] = powerspec(matrix,doPlot,corners,[tmpvar],[dofft]);
% Input:
%    matrix  = input image
%    doPlot  = plot results ('y' or 'n') [default = 'n']
%    corners = include corners ('y' or 'n') [default = 'n']
%    tmpvar  = variable to speed up next analysis
%    dofft   = do fft before plotting powerspectrum [default = 'y']

% ---------------------------------
% Copyright (c) 2004 Serge Dumoulin
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License version 2
% as published by the Free Software Foundation.
% ---------------------------------

% input args
if ~exist('matrix','var') || isempty(matrix),
    help(mfilename);
    return;
end;
if ~exist('doPlot','var') || isempty(doPlot),
    doPlot=false;
end;
if ~exist('corners','var') || isempty(corners),
    corners=false;
end;
if ~exist('dofft','var') || isempty(dofft),
    dofft = true;
end;

% if tmpvar is not given recompute some parameters otherwise load them
if ~exist('tmpvar','var'),
    recompFreq = true;

    % check
    [len,wid]=size(matrix);
    if len ~= wid,
        disp(sprintf('Error --- not a square matrix'));
    end;

    % corners?
    if corners,
        endspec = round(max(dist_from_center(:)));
    else,
        endspec = round(len/2);
    end;

else,
    recompFreq = false;
    avgmatrix  = tmpvar{1};
    keep       = tmpvar{2};
    endspec    = tmpvar{3};
end;


% get fft
if dofft,
    matrix = matrix-mean(matrix(:));
    matrix = abs(fftshift(fft2(matrix)));
    % convert to rms units
    matrix = sqrt(2).*matrix./numel(matrix);
    matrix = matrix(:);
end;


% get rel freq
if recompFreq,
    coords  = ceil(-len/2:wid/2-1);
    [xx,yy] = meshgrid(coords);
    dist_from_center = round(sqrt((xx.*xx)+(yy.*yy)));

    dist_from_center  = dist_from_center(:);
    keep              = dist_from_center<=endspec;
    dist_from_center  = dist_from_center(keep);

    % compute average matrix
    warning('off','MATLAB:divideByZero');
    freq = 0:endspec;
    avgmatrix = sparse(numel(dist_from_center),numel(freq));
    for n=1:numel(freq), %vectorize?
        avgmatrix(:,n) = dist_from_center==freq(n);
        avgmatrix(:,n) = avgmatrix(:,n)./sum(avgmatrix(:,n));
    end
    warning('on','MATLAB:divideByZero');

    if nargout > 1,
        tmpvar{1} = avgmatrix;
        tmpvar{2} = keep;
        tmpvar{3} = endspec;
    end;
end;

% limit matrix to the data we are using
matrix = matrix(keep)';

% average powerspec across orientations
powerspectrum = matrix*avgmatrix;
powerspectrum(isnan(powerspectrum))=0;
powerspectrum = powerspectrum(:);

% plot?
if doPlot,
    figure;
    ToPlot=abs(powerspectrum);
    ToPlot=ToPlot(2:length(ToPlot));
    loglog(ToPlot);axis('square');
end;

return;
