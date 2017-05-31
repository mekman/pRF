function dispPath = getDisplayPath()
% getDisplayPath - path where screen information resides
%
% dispPath = getDisplayPath
%

dispPath = [fileparts(fileparts(fileparts(which(mfilename)))) filesep 'Displays'];

if ~exist(dispPath,'dir')
    warning(sprintf('[%s]: display path does not exist (%s)\n',mfilename,dispPath));
end

return;