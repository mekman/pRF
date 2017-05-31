
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


%duration of first and second image presentations, and inter-stimulus interval. Expressed in units of 1/16th of a second    
firstOnPeriod= 4;
secondOnPeriod= 4;
ISIPeriod=8;


% load matrix or make it
if ~isempty(params.loadMatrix),
    % we should really put some checks that the matrix loaded is
    % appropriate etc.
    load(params.loadMatrix);
    disp(sprintf('[%s]:loading images from %s.',mfilename,params.loadMatrix));
%    disp(sprintf('[%s]:size stimulus: %dx%d pixels.',mfilename,n,m));
else
    outerRad = params.radius;

    %%% Set check colormap indices %%%
    %bk = findName(params.display.reservedColor,'background');
    %minCmapVal = max([params.display.reservedColor(:).fbVal])+1;
    %maxCmapVal = params.display.numColors-1;
    bk = params.display.backColorIndex;

    %%% Initialize image template %%%
    m=angle2pix(params.display,2*outerRad); 
    n=angle2pix(params.display,2*outerRad);
    
    % should really do something more intelligent, like outerRad-fix
    [x,y]=meshgrid(linspace(-outerRad,outerRad,n),linspace(outerRad,-outerRad,m));


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

    
    % Initialize empty image structure, to be filled later

    images=zeros(m,n,1,'uint8');

end;
%save 
%return;


% Now we compile the actual sequence

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

sequence = repmat(sequence(:),params.ncycles,1);


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
    firstOnB=firstOn;
    firstOnB(:)=2;

    secondOn=zeros(secondOnPeriod, 1);
    secondOnB=secondOn;
    secondOnB(:)=3;
    
    ISI=ones(ISIPeriod, 1);
    off=ones(offPeriod, 1);
    
    oneCycle=[firstOn; ISI; secondOn;off];
    oneCycleB=[firstOnB; ISI; secondOnB; off];
    
    seq2=repmat(oneCycle, params.ncycles,1);
    sequence=repmat(oneCycleB, params.ncycles, 1);
    
    seq2=logical(seq2);
    
    
    images=zeros(m,n,1,'uint8');
    if isempty(params.loadMatrix),
        sequence(seq2)    = 1;
        images(:,:,1)   = uint8(ones(size(images,1),size(images,2)).*bk);
    end;
    clear seq2;    
    
    im=loadProcessedImages(IMAGEPATH, length(images(:,1,1)));
    if length(im)>params.insertBlanks.freq
        im=im(1:params.insertBlanks.freq);
    end
    
    
     for cyclecounter=1:length(im)
        disp(sprintf('[%s]:image %.0f is %s', mfilename, cyclecounter, im(cyclecounter).filename));
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

stimulus.pictures=im;

% save matrix if requested
if ~isempty(params.saveMatrix),
    save(params.saveMatrix,'images');
end;

