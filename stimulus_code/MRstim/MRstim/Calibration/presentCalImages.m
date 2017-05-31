function presentCalImages
% presentCalImages - present mean-luminance images for monitior calibration

% 11/2005 SOD: wrote it.


% default parameters
params.duration     = 2; % seconds
params.nimages      = 9;
params.nindex       = fliplr(round(linspace(0,255,params.nimages)));%[0:256/(params.nimages-1)]*(params.nimages-1);
params.quitProgKey  = KbName('q'); 
% empty is no calibration
params.calibration  = [];
% empty is no timing but sunc to button press
params.timing       = [];%[1:length(params.nindex)*3].*params.duration;

if ~isempty(params.calibration),
    params.display = loadDisplayParams('displayName',params.calibration);
    disp(sprintf('[%s]:loading calibration from: %s.',mfilename,params.calibration));
else,
    params.display.screenNumber   = max(Screen('screens'));
    [width, height]=Screen('WindowSize',params.display.screenNumber);
    params.display.numPixels  = [width height];
    params.display.dimensions = [24.6 18.3];
    params.display.pixelSize  = min(params.display.dimensions./params.display.numPixels);
    params.display.distance   = 40;
    params.display.frameRate  = 60;
    params.display.cmapDepth  =  8;
    params.display.gammaTable = [0:255]'./255*[1 1 1];
    params.display.gamma      = params.display.gammaTable;
    params.display.backColorRgb   = [128 128 128 255];
    params.display.textColorRgb   = [255 255 255 255];
    params.display.backColorIndex = 128;
    params.display.stimRgbRange   = [0 255];
    params.display.bitsPerPixel   = 32;
    disp(sprintf('[%s]:no calibration.',mfilename));    
end;
params.display.quitProgKey = params.quitProgKey;
params.display.numColors   = size(params.display.gamma,1);

params.display.devices     = getDevices;



corder = [eye(3); ones(1,3)];
allcolors = [];
for n=1:numel(params.nindex),
    allcolors = [allcolors;corder.*params.nindex(n)];
end;
allcolors./255

try,
    % open screen
    calGamma = params.display.gamma(round(linspace(1,size(params.display.gamma,1), params.display.numColors)),:);
    params.display.oldGamma = Screen('LoadNormalizedGammaTable', params.display.screenNumber, calGamma);

    %Open the screen and save the window pointer and rect
    numBuffers = 2;
    [params.display.windowPtr,params.display.rect] = ...
        Screen('OpenWindow',params.display.screenNumber,...
               params.display.backColorRgb, [], params.display.bitsPerPixel, numBuffers);

    % Give the display a moment to recover from the change of display mode when
    % opening a window. It takes some monitors and LCD scan converters a few seconds to resync.
    WaitSecs(2);

                   
    Screen('FillRect', params.display.windowPtr, params.display.backColorRgb);
    Screen('Flip', params.display.windowPtr);
    HideCursor;

    %params.display                = openScreen(params.display);
    
    % ready 
    sound(sin(1:100));
    disp(sprintf('[%s]:Press any key to begin',mfilename));
    pause;
    KbCheck;
    
    % go
    t0 = GetSecs;
    quitProg=0;
    for n=1:size(allcolors,1),
        color = allcolors(n,:);
        disp(sprintf('[%d]:\t[%d %d %d]\t=[%f %f %f]',n,color,color./255));
        sound(sin(1:100));
        Screen('FillRect', params.display.windowPtr, color);
        %            Screen('DrawDots', params.display.windowPtr,round(params.display.numPixels'./2),min(params.display.numPixels), color);
        Screen('Flip', params.display.windowPtr);
        if ~isempty(params.timing),
            [params.timing(ind1) color]
            while getSecs-t0 < params.timing(ind1),
                [exKeyIsDown,exSecs,exKeyCode] = KbCheck(params.display.devices.keyInputInternal);
                if(exKeyIsDown),
                    if(exKeyCode(params.quitProgKey)),
                        quitProg = 1;
                        break; % out of while loop
                    end;
                end;
            end;
            if quitProg,break;end;
        else,
            pause;
        end;
    end;
    
    % Close the one on-screen and many off-screen windows
    closeScreen(params.display);

catch,
    % clean up if error occurred
    Screen('CloseAll');
    setGamma(0);
    ShowCursor;
    rethrow(lasterror);
end;

        
