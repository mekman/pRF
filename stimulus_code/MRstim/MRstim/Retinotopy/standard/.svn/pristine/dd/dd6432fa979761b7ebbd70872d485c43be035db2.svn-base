
function stimulus = makeRetinotopyStimulusImages(params)
% makeRetinotopyStimulus - make various retinotopy stimuli
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

IMAGEPATH='/Applications/Psychtoolbox/MRstim/Display images fMRI/Candidate images/Processed';


% load matrix or make it
if ~isempty(params.loadMatrix),
    % we should really put some checks that the matrix loaded is
    % appropriate etc.
    load(params.loadMatrix);
    disp(sprintf('[%s]:loading images from %s.',mfilename,params.loadMatrix));
%    disp(sprintf('[%s]:size stimulus: %dx%d pixels.',mfilename,n,m));
else,
    outerRad = params.radius;
    innerRad = params.innerRad;
    wedgeWidth = params.wedgeWidth;
    ringWidth = params.ringWidth;

    numImages = params.numImages;
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
%     if isfield(params.display, 'Rect')
%         y=y+pix2angle(params.display, params.display.Rect(2)/2);
%     end
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
    
    mask = makecircle(m*24/28,m);
    
    % r = eccentricity; theta = polar angle
    r = sqrt (x.^2  + y.^2);
    theta = atan2 (y, x);					% atan2 returns values between -pi and pi
    theta(theta<0) = theta(theta<0)+2*pi;	% correct range to be between 0 and 2*pi

    % Calculate checkerboard.
    % Wedges alternating between -1 and 1 within stimulus window.
    % The computational contortions are to avoid sign=0 for sin zero-crossings
    wedges = sign(2*round((sin(theta*numSubWedges*(2*pi/wedgeWidth))+1)/2)-1);
    posWedges = find(wedges==1);
    negWedges = find(wedges==-1);

    rings = wedges.*0;
    rings = sign(2*round((sin(r*numSubRings*(2*pi/ringWidth))+1)/2)-1);

    checks   = zeros(size(rings,1),size(rings,2),params.temporal.motionSteps);
    for ii=1:numMotSteps,
        tmprings1 = sign(2*round((sin(r*numSubRings*(2*pi/ringWidth)+(ii-1)/numMotSteps*2*pi)+1)/2)-1);
        tmprings2 = sign(2*round((sin(r*numSubRings*(2*pi/ringWidth)-(ii-1)/numMotSteps*2*pi)+1)/2)-1);
        rings(posWedges)=tmprings1(posWedges);
        rings(negWedges)=tmprings2(negWedges);

        checks(:,:,ii)=minCmapVal+ceil((maxCmapVal-minCmapVal) * (wedges.*rings+1)./2);
    end;
    
    % Loop that creates the final images
    fprintf('[%s]:Creating %d images:',mfilename,numImages);
    images=zeros(m,n,numImages*params.temporal.motionSteps,'uint8');
    for imgNum=1:params.numImages

        switch params.type;
            case 'wedge',
                switch(lower(params.display.fixType))
                    case {'left disk','right disk'},
                        loAngle = 2*pi*((imgNum-1)/(numImages*2));
                        hiAngle = loAngle + wedgeWidth;
                        % for second wedge
                        loAngle2 = loAngle + pi;
                        hiAngle2 = loAngle2 + wedgeWidth;
                    otherwise
                        loAngle = 2*pi*((imgNum-1)/numImages);
                        hiAngle = loAngle + wedgeWidth;
                end;
                loEcc = innerRad;
                hiEcc = outerRad;
            case 'ring',
                loAngle = 0;
                hiAngle = 2*pi;
                loEcc = outerRad * (imgNum-1)/numImages;
                hiEcc = loEcc+ringWidth;
            case 'center-surround',
                loAngle = 0;
                hiAngle = 2*pi;
                if mod(imgNum,2)
                    loEcc = params.centerInnerRad;
                    hiEcc = params.centerOuterRad;
                else
                    loEcc = params.surroundInnerRad;
                    hiEcc = params.surroundOuterRad;
                end
                switch(lower(params.display.fixType))
                    case {'left disk','right disk'},
                        hiEcc = hiEcc.*2;
                end;
            otherwise,
                error('Unknown stimulus type!');

        end
        % This isn't as bad as it looks
        % Can fiddle with this to clip the edges of an expanding ring - want the ring to completely 
        % disappear from view before it re-appears again in the middle.

        % Can we do this just be removing the second | from the window expression? so...
        window = ( ((theta>=loAngle & theta<hiAngle) | ...
            (hiAngle>2*pi & theta<mod(hiAngle,2*pi))) & ...
            ((r>=loEcc & r<=hiEcc)) & ...
            r<outerRad) ;
    
%         if isfield(params.display, 'Rect')
%             windowRadius = pix2angle(params.display, (params.display.Rect(4)-params.display.Rect(2)))/2;
%             window=window & r<windowRadius;
%         end
        
        % if we have a fixation to the side double the wedges
        switch(lower(params.display.fixType))
            case {'left disk','right disk'},
                switch params.type;
                    case 'wedge',
                        window = window +...
                            ( ((theta>=loAngle2 & theta<hiAngle2) | ...
                            (hiAngle2>2*pi & theta<mod(hiAngle2,2*pi))) & ...
                            ((r>=loEcc & r<=hiEcc)) & ...
                            r<outerRad);
                        window = window > .5;
                        % hack!
                end;
                window = window .* mask;
                window = window > .5;
        end;

        % yet another loop to be able to move the checks...
        for ii=1:numMotSteps, 
            img = bk*ones(m,n);
            tmpvar = checks(:,:,ii);
            img(window) = tmpvar(window);	
            images(:,:,imgNum*numMotSteps-numMotSteps+ii) = uint8(img);
        end;
        fprintf('.');drawnow;
    end
    fprintf('Done.\n');
end;
%save 
%return;


% Now we compile the actual sequence
seq = [];
curCmapFrame = 1;
if params.seqDirection==0
	imSeq = [1:params.numImages];
else
	imSeq = [params.numImages:-1:1];
end

duration.stimframe          = 1./params.temporal.frequency./params.temporal.motionSteps;
duration.scan.seconds       = params.ncycles*params.period;
duration.scan.stimframes    = params.ncycles*params.period./duration.stimframe;
duration.cycle.seconds      = params.period;
duration.cycle.stimframes   = params.period./duration.stimframe;
duration.prescan.seconds    = params.prescanDuration;
duration.prescan.stimframes = params.prescanDuration./duration.stimframe;

% make stimulus sequence
% main wedges/rings
sequence = ...
    ones(1,duration.cycle.stimframes./params.numImages)'*...
    [1:params.temporal.motionSteps:params.temporal.motionSteps*params.numImages];

if params.insertBlanks.do,
    if params.insertBlanks.phaseLock, % keep 1 cycle and repeat
        completeCycle = sequence(:);
        sequence = [completeCycle;...
                    completeCycle(1:round(end/2));...
                    completeCycle;...
                    completeCycle(1:round(end/2));...
                    completeCycle;...
                    completeCycle(1:round(end/2));...
                    completeCycle;...
                    completeCycle(1:round(end/2))];
    else,
        sequence = repmat(sequence(:),params.ncycles,1);
    end;
else,
    sequence = repmat(sequence(:),params.ncycles,1);
end;


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

% fixation dot sequence
fixSeq = ones(nn,1)*round(rand(1,ceil(length(sequence)/nn)));
fixSeq = fixSeq(:)+1;
fixSeq = fixSeq(1:length(sequence));
% force binary
fixSeq(fixSeq>2)=2; 
fixSeq(fixSeq<1)=1;

% direction
if params.seqDirection~=0
	sequence = flipud(sequence);
end

% insert blanks
if params.insertBlanks.do,
    seq2      = zeros(size(sequence));
    oneCyclePeriod  = length(seq2)/params.insertBlanks.freq;
    
    firstOnPeriod= 4;
    secondOnPeriod= 4;
    ISIPeriod=8;
    offPeriod=oneCyclePeriod-ISIPeriod-firstOnPeriod-secondOnPeriod;
    
    firstOn=zeros(firstOnPeriod, 1);
    firstOn(:)=1;
    
    secondOn=zeros(secondOnPeriod, 1);
    secondOn(:)=2;
    
    ISI=zeros(ISIPeriod, 1);
    off=zeros(offPeriod, 1);
    
    oneCycle=[firstOn; ISI; secondOn;off];

    imageslength=size(images, 3)+1;
    
    seq2=repmat(oneCycle, params.insertBlanks.freq,1);
    
    im=loadProcessedImages(IMAGEPATH, length(images(:,1,1)));
    
    cyclecounter=1;
    disp(sprintf('[%s]:image %.0f is %s', mfilename, cyclecounter, im(cyclecounter).filename));
    for i=1:imageslength
        if i>1
            if seq2(i)==1 && seq2(i-1)==0
                randomizer=round(rand(1));
                cyclecounter=cyclecounter+1;
                if cyclecounter>length(im)
                    cyclecounter=1;
                end
                disp(sprintf('[%s]:image %.0f is %s', mfilename, cyclecounter, im(cyclecounter).filename));
            end
        else
            randomizer=round(rand(1));
        end
        if seq2(i)==1
            if randomizer==1
                images(:,:,i)=im(cyclecounter).circleimage1;
            else
                images(:,:,i)=im(cyclecounter).circleimage2;
            end
        elseif seq2(i)==2
            if randomizer==0
                images(:,:,i)=im(cyclecounter).circleimage1;
            else
                images(:,:,i)=im(cyclecounter).circleimage2;
            end
        elseif seq2(i)==0
            images(:,:,i)=uint8(ones(size(images,1),size(images,2)).*bk);
        end
    end

    disp(sprintf('[%s]:Stimulus on for %.1f and off for %.1f seconds.',...
        mfilename,(firstOnPeriod+secondOnPeriod)*duration.stimframe,(ISIPeriod+offPeriod)*duration.stimframe));
end;
        
% Insert the preappend images by copying some images from the
% end of the seq and tacking them on at the beginning
sequence = [sequence(length(sequence)+1-duration.prescan.stimframes:end); sequence];
timing   = [0:length(sequence)-1]'.*duration.stimframe;
cmap     = params.display.gammaTable;
fixSeq   = [fixSeq(length(fixSeq)+1-duration.prescan.stimframes:end); fixSeq];




% make stimulus structure for output
stimulus = createStimulusStruct(images,cmap,sequence,[],timing,fixSeq);

% save matrix if requested
if ~isempty(params.saveMatrix),
    save(params.saveMatrix,'images');
end;

