function displayID = openScreen(displayID, hideCursorFlag)
% openScreen - open screen for psychtoolbox
%
% Usage: displayID = openScreen(displayID, [hideCursorFlag=true])
%
% openScreen takes a displayID structure (e.g., from running
% loadDisplayParams), and does the following:
% 1. opens a PTB window with a background color (defaults to 1/2 of
% maxRgbValue) (using Screen in PTB)
% 2. draws a fixation dot (using drawFixation in exptTools)
% 3. hides the cursor
% 4. store the original gamma table in the displayID structure
%
% After you are done with the opened PTB window, use closeScreen to revert
% back to the original state.
%
% History:
% ##/##/## rfd & sod wrote it.
% 04/12/06 shc (shcheung@stanford.edu) cleaned it and added the help
% comments.

if(~exist('hideCursorFlag','var')||isempty(hideCursorFlag))
    hideCursorFlag = true;
end

if(~isfield(displayID,'screenNumber'))
    displayID.screenNumber = 0;
    disp('Using default screen number of 0');
end

if(~isfield(displayID,'frameRate'))
    displayID.frameRate = 60;
    disp('Using default frame rate of 60 Hz');
end

if(~isfield(displayID,'bitsPerPixel'))
    displayID.bitsPerPixel = 8;
    disp('Using default pixel depth of 8 bits');
end

if(~isfield(displayID,'pixelSize'))
    if(isfield(displayID,'numPixels') && isfield(displayID,'dimensions') && ...
            ~isempty(displayID.numPixels) && ~isempty(displayID.dimensions) && ...
            length(displayID.numPixels) == length(displayID.dimensions))
        displayID.pixelSize = mean(displayID.dimensions./displayID.numPixels);
        disp('Using number of pixels and dimension information to calculate pixel size');
    else
        displayID.pixelSize = 0.0691;
        disp('Using default pixel size of 0.0691 cm');
    end
end

if(~isfield(displayID,'gammaTable'))
    displayID.gammaTable = linspace(0,2^displayID.cmapDepth-1,2^displayID.cmapDepth)';
    disp('Using default linear gamma table');
end

if(~isfield(displayID,'distance'));
    displayID.distance = 50;
    disp('Using default display distance of 50 cm');
end

if(~isfield(displayID,'backColorRgb'))
    displayID.backColorRgb = [repmat(round(displayID.maxRgbValue/2),1,3) displayID.maxRgbValue];
    disp(['Setting backColorRgb to ',num2str(displayID.backColorRgb),'.']);
end

% Skip the annoying blue flickering warning
Screen('Preference','SkipSyncTests',0);

displayID.oldGamma = Screen('ReadNormalizedGammaTable', displayID.screenNumber);
try
    Screen('LoadNormalizedGammaTable', displayID.screenNumber,displayID.gamma);
catch
    % displayID.gamma may be 10bit in that case reduce to 8 spanning entire
    % range
    fprintf(1,'[%s]:error load 10bit lookup table, reverting to 8bit.',mfilename);
    putgamma = displayID.gamma(round(linspace(1,size(displayID.gamma,1),256)),:);
    Screen('LoadNormalizedGammaTable', displayID.screenNumber,putgamma);
end;

% Open the screen and save the window pointer and rect
numBuffers = 2;
[displayID.windowPtr,displayID.rect] = Screen('OpenWindow',displayID.screenNumber,...
    displayID.backColorRgb, [], displayID.bitsPerPixel, numBuffers);

Screen('FillRect', displayID.windowPtr, displayID.backColorRgb);
drawFixation(displayID);
Screen('Flip', displayID.windowPtr, 3);
if(hideCursorFlag)
    HideCursor;
end
return;
