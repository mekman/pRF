function estimateMidGray
% estimateMidGray - estimate mid gray DAC value

% 2008/06 SOD: wrote it

oldLevel = Screen('Preference', 'Verbosity', 1);

% default parameters
params.duration     = 2;
params.nimages      = 5;
params.nindex       = fliplr(round(linspace(0,255,5)));%[0:256/(params.nimages-1)]*(params.nimages-1);
params.quitProgKey  = KbName('q'); 
% empty is no calibration
params.calibration  = '3T_projector_UMC_800x600';,%[];
% empty is no timing but sunc to button press
params.timing       = [];%[1:length(params.nindex)*3].*params.duration;

if ~isempty(params.calibration),
    params.display = loadDisplayParams('displayName',params.calibration);
    fprintf(1,'[%s]:loading calibration from: %s.\n',mfilename,params.calibration);
else,
    params.display.screenNumber   = max(Screen('screens'));
    [width, height]=Screen('WindowSize',params.display.screenNumber);
    params.display.numPixels  = [width height];
    params.display.cmapDepth  =  8;
    params.display.gammaTable = [0:255]'./255*[1 1 1];
    params.display.gamma      = params.display.gammaTable;
    params.display.backColorRgb   = [128 128 128 255];
    params.display.textColorRgb   = [255 255 255 255];
    params.display.backColorIndex = 128;
    params.display.stimRgbRange   = [0 255];
    params.display.bitsPerPixel   = 32;
    params.display.fixType        = 'disk';
    params.display.fixColorRgb    = [255 0 0];
    params.display.fixX           = width./2;
    params.display.fixY           = height./2;
    params.display.fixSizePixels  = 4;
    
    disp(sprintf('[%s]:no calibration.',mfilename));    
end;
params.display.quitProgKey = params.quitProgKey;
params.display.numColors   = size(params.display.gamma,1);

params.display.devices     = getDevices(false);

% make calibration image
LineWidth = 2;
calline1  = rem(round((1:params.display.numPixels(1))/LineWidth),2).*2-1;
calline2  = rem(round((1:params.display.numPixels(1))/(params.display.numPixels(1)/10)),2);
calline   = calline1.*calline2;
calline   = round((calline/2 + .5).*255);
stimulus.images  = ones(params.display.numPixels(2),1)*calline;


key.up   = KbName('u');
key.down = KbName('d');
key.end  = KbName('q');

quitProg = false;
gamma = params.display.gamma.*255;
try,
    % open screen
    params.display = openScreen(params.display);
    
    % Give the display a moment to recover from the change of display mode when
    % opening a window. It takes some monitors and LCD scan converters a few seconds to resync.
    WaitSecs(2);

    Screen('FillRect', params.display.windowPtr, params.display.backColorRgb);

    % load stimulus
    stimulus = createTextures(params.display, stimulus);
    Screen('DrawTexture', params.display.windowPtr, stimulus.textures(1), stimulus.srcRect, stimulus.destRect);


    Screen('Flip', params.display.windowPtr);
    HideCursor;

    
    % ready 
    sound(sin(1:100));
    KbCheck;
    
    % go
    Screen('DrawTexture', params.display.windowPtr, stimulus.textures(1), stimulus.srcRect, stimulus.destRect);
    fprintf(1,'[%s]:Initial mid gray DAC value is [%d %d %d]\n',mfilename,gamma(129,:));
    while(1),
        % scan the keyboard for experimentor input
        [exKeyIsDown,exSecs,exKeyCode] = KbCheck(params.display.devices.keyInputInternal);
        if(exKeyIsDown)
            if(exKeyCode(key.end)),
                quitProg = true;
                fprintf(1,'[%s]:Mid gray DAC values is estimated as [%d %d %d].\n\n',mfilename,gamma(129,:));
                break; % out of while loop
            end;
            if(exKeyCode(key.up)),
                gamma(129,:) = gamma(129,:)+1;
            end;
            if(exKeyCode(key.down)),
                gamma(129,:) = gamma(129,:)-1;
            end;            
            gamma=min(gamma,255);
            gamma=max(gamma,0);
            Screen('LoadNormalizedGammaTable', params.display.screenNumber,gamma./255);
        end;
        % release CPU time
        WaitSecs(0.01);
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

Screen('Preference', 'Verbosity', oldLevel);
