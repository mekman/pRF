function [stimulus otherSequence] = makeRetinotopyStimulus_bars8Pass(params)
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


blanksAfterPass=2;  %1 (True) puts blanks after bars have passed, every pass. 0 (false) puts blanks within the bar pass. 2 puts blanks after only horizontal bar passes, which is what we do
blankTime=10.5;       %Blank time (in seconds) will be rounded up to the nearest TR
barMotionReversalTime=2; %Checkerboard motion direction is parallel to bar orientation, and reverses at the most this frequently (randomized)
fixationChangeTime=2;    %Fixation dot changes color (green-red or red-green) at the most this frequently (randomized)


% various time measurements:
duration.stimframe          = 1./params.temporal.frequency./params.temporal.motionSteps;
duration.scan.seconds       = params.ncycles*params.period;
duration.scan.stimframes    = params.ncycles*params.period./duration.stimframe;
duration.cycle.seconds      = params.period;
duration.cycle.stimframes   = params.period./duration.stimframe;
duration.prescan.seconds    = params.prescanDuration;
duration.prescan.stimframes = params.prescanDuration./duration.stimframe;

if blanksAfterPass==2;
    totalBlankTime   = 4.*ceil(blankTime./params.tr).*params.tr;
    barPassTime=duration.cycle.seconds-totalBlankTime;
    numBarImages=params.numImages-(totalBlankTime./params.tr);
    barPassFrames=barPassTime/duration.stimframe;
elseif blanksAfterPass==1;
    totalBlankTime   = 8.*ceil(blankTime./params.tr).*params.tr;
    barPassTime=duration.cycle.seconds-totalBlankTime;
    numBarImages=params.numImages-(totalBlankTime./params.tr);
    barPassFrames=barPassTime/duration.stimframe;
else
    barPassTime=duration.cycle.seconds;
    numBarImages=params.numImages;
    barPassFrames=duration.cycle.stimframes;
end


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

    halfNumImages = numBarImages./2;
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
    orientations = orientations([3 2 5 4 7 6 1 8]);  %[6 3 8 5 2 7 4 1]);
    remake_xy    = zeros(1,numBarImages)-1;
    remake_xy(1:length(remake_xy)/length(orientations):length(remake_xy)) = orientations;
    original_x   = x;
    original_y   = y;
    % step size of the bar
    step_nx      = barPassTime./params.tr/8;
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

                tmpvar = zeros(m,n);
                tmpvar(window) = 1;
                tmpvar = repmat(tmpvar,[1 1 numMotSteps]);
                window = tmpvar == 1;
                img         = bk*ones(size(checks));
                img(window) = checks(window);
                images(:,:,(imgNum-1).*numMotSteps+1:imgNum.*numMotSteps) = uint8(img);


        fprintf('.');drawnow;
    end
    fprintf('Done.\n');
end;



% Now we compile the actual sequence
% make stimulus sequence, make half and then add the rest as a flipped
% version of the first half
sequence = ...
    ones(barPassFrames./2./halfNumImages,1)*...
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
%nn=30; % this should be a less random choice, ie in seconds
nn = 1./duration.stimframe*barMotionReversalTime; % on average every 4 seconds [max response time = 3 seconds]
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

% insert blanks
if params.insertBlanks.do,
    if blanksAfterPass==2
        blankCycleFrames=length(sequence)/(params.insertBlanks.freq*2);
        offTime   = ceil(blankTime./params.tr).*params.tr; % make sure it's a multiple of the tr
        offPeriod = ceil(offTime./duration.stimframe);
        
        images(:,:,size(images,3)+1)   = uint8(ones(size(images,1),size(images,2)).*bk);
        blankInsert=ones(offPeriod,1).*size(images,3);
        sequence=cat(1,sequence(1:blankCycleFrames), blankInsert, sequence(blankCycleFrames+1:3*blankCycleFrames), blankInsert, sequence(blankCycleFrames*3+1:5*blankCycleFrames),blankInsert, sequence(blankCycleFrames*5+1:7*blankCycleFrames), blankInsert, sequence(blankCycleFrames*7+1:end));

    elseif blanksAfterPass==1
        blankCycleFrames=length(sequence)/length(Orientations);
        offTime   = ceil(blankTime./params.tr).*params.tr; % make sure it's a multiple of the tr
        offPeriod = ceil(offTime./duration.stimframe);
        
        images(:,:,size(images,3)+1)   = uint8(ones(size(images,1),size(images,2)).*bk);
        blankInsert=ones(offPeriod,1).*size(images,3);
        newSequence=[];
        for passCounter=1:length(Orientations)
            newSequence=cat(1, newSequence, sequence((passCounter-1)*blankCycleFrames+1:passCounter*blankCycleFrames), blankInsert);
        end
        sequence=newSequence;

    else
        seq2      = zeros(size(sequence));
        oneCycle  = length(seq2)/params.insertBlanks.freq;
        offTime   = ceil(blankTime./params.tr).*params.tr; % make sure it's a multiple of the tr
        offPeriod = ceil(offTime./duration.stimframe);
        onPeriod  = oneCycle-offPeriod;
        oneBlankCycle = [false(onPeriod,1); true(offPeriod,1)]; % blanks at the end
        oneBlankCycle = [oneBlankCycle(oneCycle/2+1:end); oneBlankCycle(1:oneCycle/2)]; % blanks half way
        seq2      = repmat(oneBlankCycle,params.insertBlanks.freq,1);
        add       = size(images,3)+1;
        if isempty(params.loadMatrix),
            sequence(seq2)    = add;
            images(:,:,add)   = uint8(ones(size(images,1),size(images,2)).*bk);
        end;
        clear seq2;
        disp(sprintf('[%s]:Stimulus on for %.1f and off for %.1f seconds.',...
            mfilename,onPeriod*duration.stimframe,offPeriod*duration.stimframe));
    end
end;
        

% fixation dot sequence
% change on the fastest every 6 seconds
minsec = fixationChangeTime./duration.stimframe;
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

otherSequence = [];

% make stimulus structure for output
stimulus = createStimulusStruct(images,cmap,sequence,[],timing,fixSeq);

% save matrix if requested
if ~isempty(params.saveMatrix),
    save(params.saveMatrix,'images');
end;

