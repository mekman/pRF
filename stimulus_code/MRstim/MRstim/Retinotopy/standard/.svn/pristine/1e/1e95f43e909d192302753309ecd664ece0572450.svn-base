function [stimulus otherSequence] = makeRetinotopyStimulus_bars(params)
% makeRetinotopyStimulus - make various retinotopy stimuli
%
% stimulus = makeRetinotopyStimulus_bars(params)
%
% Matlab code to generate various retinotopy stimuli
% Generates one full cycle, as well as the sequence for the entire scan.
%
% 99.09.15 RFD: I fixed the sequence generation algorithm so that
%   timing is now frame-accurate.  The algorithm now keeps track
%   of timing error that accumulates due to rounding to the nearest
%   frame and corrects for that error when it gets to be more than 
%   half a frame.  
%   The algorithm also randomely reverses the drift direction, rather
%   than reversing every half-an image duration.
% 2005.06.15 SOD: changed for OSX - stimulus presentation will now be 
%                 time-based rather than frame based. Because of bugs
%                 with framerate estimations.


% various time measurements:
duration.stimframe          = 1./params.temporal.frequency./params.temporal.motionSteps;
duration.scan.seconds       = params.ncycles*params.period;
duration.scan.stimframes    = params.ncycles*params.period./duration.stimframe;
duration.cycle.seconds      = params.period;
duration.cycle.stimframes   = params.period./duration.stimframe;
duration.prescan.seconds    = params.prescanDuration;
duration.prescan.stimframes = params.prescanDuration./duration.stimframe;


% load matrix or make it
if ~isempty(params.loadMatrix),
    % we should really put some checks that the matrix loaded is
    % appropriate etc.
    load(params.loadMatrix);
    halfNumImages = params.numImages./2;
    disp(sprintf('[%s]:loading images from %s.',mfilename,params.loadMatrix));
%    disp(sprintf('[%s]:size stimulus: %dx%d pixels.',mfilename,n,m));
else
    outerRad = params.radius;
    innerRad = params.innerRad;
    wedgeWidth = params.wedgeWidth;
    ringWidth = params.ringWidth;

    halfNumImages = params.numImages./2;
    numMotSteps = params.temporal.motionSteps;
    numSubRings = params.numSubRings;
    numSubWedges = params.numSubWedges;

    %%% Set check colormap indices %%%
    %bk = findName(params.display.reservedColor,'background');
    %minCmapVal = max([params.display.reservedColor(:).fbVal])+1;
    %maxCmapVal = params.display.numColors-1;
    bk = params.display.backColorIndex;
    
    minCmapVal = min([params.display.stimRgbRange]);
    maxCmapVal = max([params.display.stimRgbRange]);


    %%% Initialize image template %%%
    m=angle2pix(params.display,2*outerRad); 
    n=angle2pix(params.display,2*outerRad);
    
    % should really do something more intelligent, like outerRad-fix
    switch(lower(params.display.fixType))
        case 'left disk',
            [x,y]=meshgrid(linspace( 0,outerRad*2,n),linspace(outerRad,-outerRad,m));
            outerRad = outerRad.*2;
        case 'right disk',
            [x,y]=meshgrid(linspace(-outerRad*2,0,n),linspace(outerRad,-outerRad,m));
            outerRad = outerRad.*2;
        otherwise,
            [x,y]=meshgrid(linspace(-outerRad,outerRad,n),linspace(outerRad,-outerRad,m));
    end;
    
    % here we crop the image if it is larger than the screen 
    % seems that you have to have a square matrix, bug either in my or
    % psychtoolbox' code - so we make it square
    if m>params.display.numPixels(2),
        start  = round((m-params.display.numPixels(2))/2);
        len    = params.display.numPixels(2);
        y = y(start+1:start+len, start+1:start+len);
        x = x(start+1:start+len, start+1:start+len);
        m = len;
        n = len;
    end;
    disp(sprintf('[%s]:size stimulus: %dx%d pixels.',mfilename,n,m));

    % r = eccentricity; theta = polar angle
    r = sqrt (x.^2  + y.^2);
    theta = atan2 (y, x);					% atan2 returns values between -pi and pi
    theta(theta<0) = theta(theta<0)+2*pi;	% correct range to be between 0 and 2*pi

    
    % loop over different orientations and make checkerboard
    % first define which orientations
    orientations = (0:45:360)./360*(2*pi); % degrees -> rad
    orientations = orientations([1 6 3 8 5 2 7 4]);
    remake_xy    = zeros(1,params.numImages)-1;
    remake_xy(1:length(remake_xy)/length(orientations):length(remake_xy)) = orientations;
    original_x   = x;
    original_y   = y;
    % step size of the bar
    step_nx      = duration.cycle.seconds./params.tr/8;
    step_x       = (2*outerRad) ./ step_nx;
    step_startx  = (step_nx-1)./2.*-step_x - (ringWidth./2);
    %[0:step_nx-1].*step_x+step_startx+ringWidth./2
    disp(sprintf('[%s]:stepsize: %f degrees.',mfilename,step_x));
    
    % if we create colored bars we want to make the edges soft.
    switch params.experiment,
        case {'8 bars (LMS) with blanks'}
            edgewidth = 25; % .25cpd on a 600 pixel(3deg screen)
            softmask  = makecircle(m-2*edgewidth,m,edgewidth);
        otherwise,
            softmask = ones(m);
    end;

    % Loop that creates the final images
    fprintf('[%s]:Creating %d images:',mfilename,halfNumImages);
    images=zeros(m,n,halfNumImages*params.temporal.motionSteps,'uint8');
    for imgNum=1:halfNumImages
        
        if remake_xy(imgNum) >=0,
            x = original_x .* cos(remake_xy(imgNum)) - original_y .* sin(remake_xy(imgNum));
            y = original_x .* sin(remake_xy(imgNum)) + original_y .* cos(remake_xy(imgNum));
            % Calculate checkerboard.
            % Wedges alternating between -1 and 1 within stimulus window.
            % The computational contortions are to avoid sign=0 for sin zero-crossings
            switch params.experiment
                case {'8 bars','8 bars with blanks','8 bars (slow)','8 bars with blanks (attn)'}
                    wedges    = sign(round((cos((x+step_startx)*numSubRings*(2*pi/ringWidth)))./2+.5).*2-1);
                    posWedges = find(wedges== 1);
                    negWedges = find(wedges==-1);
                    rings     = zeros(size(wedges));

                    checks    = zeros(size(rings,1),size(rings,2),params.temporal.motionSteps);
                    for ii=1:numMotSteps,
                        tmprings1 = sign(2*round((cos(y*numSubRings*(2*pi/ringWidth)+(ii-1)/numMotSteps*2*pi)+1)/2)-1);
                        tmprings2 = sign(2*round((cos(y*numSubRings*(2*pi/ringWidth)-(ii-1)/numMotSteps*2*pi)+1)/2)-1);
                        rings(posWedges) = tmprings1(posWedges);
                        rings(negWedges) = tmprings2(negWedges);

                        checks(:,:,ii)=minCmapVal+ceil((maxCmapVal-minCmapVal) * (wedges.*rings+1)./2);
                    end;
                case {'8 bars (sinewave)','8 bars (LMS)','8 bars (LMS) with blanks'}
                    
                otherwise,
                    disp(sprintf('[%s]:unknown experiment: %s.',mfilename,params.experiment));
                    return;
            end;

            % reset starting point
            loX = step_startx - step_x;
        end;


        switch params.type;
            case 'bar'
                loEcc = innerRad;
                hiEcc = outerRad;
                loX   = loX + step_x;
                hiX   = loX + ringWidth;
            otherwise,
                error('Unknown stimulus type!');

        end
        % This isn't as bad as it looks
        % Can fiddle with this to clip the edges of an expanding ring - want the ring to completely 
        % disappear from view before it re-appears again in the middle.

        % Can we do this just be removing the second | from the window
        % expression? so...
        window = ( (x>=loX & x<=hiX) & r<outerRad);
 
        % yet another loop to be able to move the checks...
        switch params.experiment
            case {'8 bars','8 bars with blanks','8 bars (slow)','8 bars with blanks (attn)'}
               
                tmpvar = zeros(m,n);
                tmpvar(window) = 1;
                tmpvar = repmat(tmpvar,[1 1 numMotSteps]);
                window = tmpvar == 1;
                img         = bk*ones(size(checks));
                img(window) = checks(window);
                images(:,:,(imgNum-1).*numMotSteps+1:imgNum.*numMotSteps) = uint8(img);
            case {'8 bars (sinewave)','8 bars (LMS)','8 bars (LMS) with blanks'}
                img1        = bk*ones(m,n);
                img2        = img1;
                tmpvar      = sin((x - loX)*numSubRings*(2*pi/ringWidth));
                img1(window) = minCmapVal+ceil((maxCmapVal-minCmapVal) * (tmpvar(window)+1)./2);
                img2(window) = minCmapVal+ceil((maxCmapVal-minCmapVal) * ((-1.*tmpvar(window))+1)./2);
                start  = imgNum*numMotSteps-numMotSteps;
%                half   = floor(numMotSteps./2);
%                images(:,:,[1:half]+start)             = repmat(uint8(img1),[1 1 half]);
%                images(:,:,[half+1:numMotSteps]+start) = repmat(uint8(img2),[1 1 numMotSteps-half]);
                
                c = sin(linspace(0,2*pi,numMotSteps+1));
                for iii = 1:numMotSteps,
                    images(:,:,iii+start) = uint8((img2-bk).*c(iii).*softmask+bk);
                end
        end

        fprintf('.');drawnow;
    end
    fprintf('Done.\n');
end;



% Now we compile the actual sequence
% make stimulus sequence, make half and then add the rest as a flipped
% version of the first half
sequence = ...
    ones(duration.cycle.stimframes./2./halfNumImages,1)*...
    [1:params.temporal.motionSteps:params.temporal.motionSteps*halfNumImages];
sequence = sequence(:);
if params.insertBlanks.do,
    if params.insertBlanks.phaseLock, % keep 1 cycle and repeat
        completeCycle = sequence;
        sequence = [completeCycle;...
                    completeCycle(1:round(end/2));...
                    completeCycle;...
                    completeCycle(1:round(end/2));...
                    completeCycle;...
                    completeCycle(1:round(end/2));...
                    completeCycle;...
                    completeCycle(1:round(end/2))];
    else,
        sequence = repmat(sequence,params.ncycles,1);
    end;
else,
    sequence = repmat(sequence,params.ncycles,1);
end;

% we make only half so we need to flip the rest
sep   = round(linspace(1,length(sequence)+1,5));
rev = []; 
for n=1:4,
    rev = [rev; flipud(sequence(sep(n):sep(n+1)-1))];
end;
sequence = [sequence; rev];



% motion frames within wedges/rings - lowpass
nn=30; % this should be a less random choice, ie in seconds
motionSeq = ones(nn,1)*round(rand(1,ceil(length(sequence)/nn)));
motionSeq = motionSeq(:)-0.5;
motionSeq = motionSeq(1:length(sequence));
motionSeq = cumsum(sign(motionSeq));

% wrap
above = find(motionSeq>params.temporal.motionSteps);
while ~isempty(above),
    motionSeq(above)=motionSeq(above)-params.temporal.motionSteps;
    above = find(motionSeq>params.temporal.motionSteps);
end;
below = find(motionSeq<1);
while ~isempty(below),
    motionSeq(below)=motionSeq(below)+params.temporal.motionSteps;
    below = find(motionSeq<1);
end;
sequence=sequence+motionSeq-1;

% direction
if params.seqDirection~=0
	sequence = flipud(sequence);
end

% insert blanks (always of for 12 seconds)
if params.insertBlanks.do,
    seq2      = zeros(size(sequence));
    oneCycle  = length(seq2)/params.insertBlanks.freq;
    offTime   = ceil(12./params.tr).*params.tr; % make sure it's a multiple of the tr
    offPeriod = ceil(offTime./duration.stimframe);
    onPeriod  = oneCycle-offPeriod;
    seq2      = repmat([zeros(onPeriod,1); ones(offPeriod,1)],params.insertBlanks.freq,1);
    add       = size(images,3)+1;
    if isempty(params.loadMatrix),
        sequence(seq2==1) = add;
        images(:,:,add)   = uint8(ones(size(images,1),size(images,2)).*bk);
    end;
    clear seq2;
    disp(sprintf('[%s]:Stimulus on for %.1f and off for %.1f seconds.',...
        mfilename,onPeriod*duration.stimframe,offPeriod*duration.stimframe));
end;
        

% fixation dot sequence
% change on the fastest every 6 seconds
minsec = 6./duration.stimframe;
fixSeq = ones(minsec,1)*round(rand(1,ceil(length(sequence)/minsec)));
fixSeq = fixSeq(:)+1;
fixSeq = fixSeq(1:length(sequence));
% force binary
fixSeq(fixSeq>2)=2; 
fixSeq(fixSeq<1)=1;


% Insert the preappend images by copying some images from the
% end of the seq and tacking them on at the beginning
sequence = [sequence(length(sequence)+1-duration.prescan.stimframes:end); sequence];
timing   = [0:length(sequence)-1]'.*duration.stimframe;
cmap     = params.display.gammaTable;
fixSeq   = [fixSeq(length(fixSeq)+1-duration.prescan.stimframes:end); fixSeq];


switch params.experiment
    case '8 bars with blanks (attn)'
        minsec = 8./duration.stimframe;
        barSeq = ones(minsec,1)*round(rand(1,ceil(length(sequence)/minsec)));
        barSeq = barSeq(:);
        barSeq = barSeq(1:length(sequence));
        % force binary
        barSeq(fixSeq>2)=2; 
        barSeq(fixSeq<1)=1;
        
        % bar flickers to lower contrast
        barSeqId = abs(barSeq(:)-circshift(barSeq(:),4))>0;
        barSeq = ones(size(barSeq));
        barSeq(barSeqId) = 2;
        barSeq=barSeq.*[0; double(diff(barSeq)~=0)]; % only keep where we need to change
        
        % mask so no change during mean luminance 
        barSeq(sequence==max(sequence))=0;
        
        % fixation flickers to lower contrast
        fixSeq = ones(minsec,1)*round(rand(1,ceil(length(sequence)/minsec)));
        fixSeq = fixSeq(:)+1;
        fixSeq = fixSeq(1:length(sequence));
        % force binary
        fixSeq(fixSeq>2)=2;
        fixSeq(fixSeq<1)=1;

        fixSeqId = abs(fixSeq(:)-circshift(fixSeq(:),8))>0;
        fixSeq = ones(size(fixSeq));
        fixSeq(fixSeqId) = 2;
        fixSeq = [fixSeq fixSeq]';
        fixSeq = fixSeq(:);
        
        % double sequence with bar contrast sequence
        sequence = [sequence -barSeq]';
        sequence = sequence(:);
        
        % put in new cmap
        cmap(:,:,2) = cmap;
        
        % images
        images(images<1) = 1; % 1-254
        images(images>254) = 254; % 1-254
        cmap(255,:,:) = cmap(256,:,:);
        cmap(255,:,2) = cmap(round(128+128.*params.attn.contrast(2,1)),:,1);
        cmap(  2,:,:) = cmap(  1,:,:);
        cmap(  2,:,2) = cmap(round(128-128.*params.attn.contrast(2,1)),:,1);
        
        % timing
        timing = [timing timing+duration.stimframe./2]';
        timing = timing(:);
        
        % output
        otherSequence = double(sequence==-2);
        
    otherwise
        otherSequence = [];
        
end



% make stimulus structure for output
stimulus = createStimulusStruct(images,cmap,sequence,[],timing,fixSeq);

% save matrix if requested
if ~isempty(params.saveMatrix),
    save(params.saveMatrix,'images');
end;

