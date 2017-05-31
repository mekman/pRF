function [response, timing, quitProg] = showScanStimulus(display,stimulus, t0)
% [response, timing] = showStimulus(display,stimulus, time0)
%
% t0 is the time the scan started and the stimulus timing should be
% relative to t0. If t0 does not exist it is created at the start of
% this program.
%  
% HISTORY:
% 2005.02.23 RFD: ported from showStimulus.
% 2005.06.15 SOD: modified for OSX. Use internal clock for timing rather
% than framesyncing because getting framerate does not always work. Using
% the internal clock will also allow some "catching up" if stimulus is
% delayed for whatever reason. Loading mex functions is slow, so this 
% should be done before callling this program.

% input checks
if nargin < 2,
	help(mfilename);
    return;
end;
if nargin < 3 || isempty(t0),
    t0 = GetSecs; % "time 0" to keep timing going
end;

% some more checks
if ~isfield(stimulus,'textures')
	% Generate textures for each image
	disp('WARNING: Creating textures before stimulus presentation.');
	disp(['         This should be done before calling ' mfilename ' for']);
	disp('         accurate timing.  See "makeTextures" for help.');
	stimulus = makeTextures(display,stimulus);
end;

% quit key
try 
    quitProgKey = display.quitProgKey;
catch
    quitProgKey = KbName('q');
end;

% some variables
nFrames = length(stimulus.seq);
HideCursor;
nGamma = size(stimulus.cmap,3);
nImages = length(stimulus.textures);
response.keyCode = zeros(length(stimulus.seq),1); % get 1 buttons max
response.secs = zeros(size(stimulus.seq));        % timing
quitProg = 0;

if isfield(stimulus, 'pictures')
    pictureCounter=0;
    nImages=nImages-length(stimulus.pictures).*2;
end

% go
disp(sprintf('[%s]:Running. Hit %s to quit.',mfilename,KbName(quitProgKey)));
for frame = 1:nFrames	
    
    %--- update display
    % If the sequence number is positive, draw the stimulus and the
    % fixation.  If the sequence number is negative, draw only the
    % fixation.

    if stimulus.seq(frame)>0
        % put in an image
		imgNum = mod(stimulus.seq(frame)-1,nImages)+1;
        if isfield(stimulus, 'pictures')
            if stimulus.seq(frame)==2
                if frame==1
                    pictureCounter=pictureCounter+1;
                    randomizer=round(rand(1));
                elseif stimulus.seq(frame-1)~=2
                    pictureCounter=pictureCounter+1;
                    randomizer=round(rand(1));                    
                end
                if nImages+pictureCounter>length(stimulus.textures)
                    pictureCounter=1;
                end
                if randomizer==0;
                    Screen('DrawTexture', display.windowPtr, stimulus.textures(nImages+pictureCounter), stimulus.srcRect, stimulus.destRect); 
                else
                    Screen('DrawTexture', display.windowPtr, stimulus.textures(nImages+pictureCounter+1), stimulus.srcRect, stimulus.destRect); 
                end
            elseif stimulus.seq(frame)==3
                if stimulus.seq(frame-1)~=3
                    pictureCounter=pictureCounter+1;
                end
                if randomizer==0;
                    Screen('DrawTexture', display.windowPtr, stimulus.textures(nImages+pictureCounter), stimulus.srcRect, stimulus.destRect); 
                else
                    Screen('DrawTexture', display.windowPtr, stimulus.textures(nImages+pictureCounter-1), stimulus.srcRect, stimulus.destRect); 
                end
            else
                Screen('DrawTexture', display.windowPtr, stimulus.textures(nImages), stimulus.srcRect, stimulus.destRect)
            end
        else            
            Screen('DrawTexture', display.windowPtr, stimulus.textures(imgNum), stimulus.srcRect, stimulus.destRect);
        end
        drawFixation(display,stimulus.fixSeq(frame));
    elseif stimulus.seq(frame)<0
        % put in a color table
		gammaNum = mod(-stimulus.seq(frame)-1,nGamma)+1;
        % The second argument is the color index.  This apparently changed
        % in recent times (07.14.2008). So, for now we set it to 1.  It may
        % be that this hsould be 
        drawFixation(display,stimulus.fixSeq(frame));
		Screen('LoadNormalizedGammaTable', display.windowPtr, stimulus.cmap(:,:,gammaNum));
    end;
    
    %--- timing
    waitTime = (GetSecs-t0)-stimulus.seqtiming(frame);
    
    %--- get inputs (subject or experimentor)
    while(waitTime<0),
        % Scan the UMC device for subject response
        [ssKeyCode,ssSecs] = deviceUMC('response',display.devices.UMCport);
        if(ssKeyCode(1)~=0)
%            kc = find(ssKeyCode);
%            response.keyCode(frame) = kc(1); % binary response for now
            response.keyCode(frame) = ssKeyCode(end); 
            response.secs(frame)    = ssSecs - t0;
        end;
        % scan the keyboard for experimentor input
        [exKeyIsDown,exSecs,exKeyCode] = KbCheck(display.devices.keyInputInternal);
        if(exKeyIsDown)
            if(exKeyCode(quitProgKey)),
                quitProg = 1;
                break; % out of while loop
            end;
        end;

        % if there is time release cpu
        if(waitTime<-0.02),
            WaitSecs(0.01);
        end;
        
        % timing
        waitTime = (GetSecs-t0)-stimulus.seqtiming(frame);
    end;
    
    %--- stop?
    if quitProg,
        disp(sprintf('[%s]:Quit signal recieved.',mfilename));
        break;
    end;

    %--- update screen
    Screen('Flip',display.windowPtr);
end;

% that's it
ShowCursor;
timing = GetSecs-t0;
disp(sprintf('[%s]:Stimulus run time: %f seconds [should be: %f].',mfilename,timing,max(stimulus.seqtiming)));

return;
