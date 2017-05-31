function [stimulus] = makeRetinotopyStimulus_barsRess(params)
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
% 2010.02.19 BMH: modified script for needs of Ress lab

%PARAMETERS FOR RESS LAB TO CHANGE AS NEEDED
switch params.experiment
    case 'Translating Bars 1', 
        barOrientations=[0 30 60 90 120 150];
    case 'Translating Bars 2', 
        barOrientations=[10 40 70 100 130 160];
    case 'Translating Bars 3',
        barOrientations=[20 50 80 110 140 170];
end

blanksAfterPass=1;  %1 (True) puts blanks after bars have passed. 0 (false) puts blanks within the bar pass.
blankTime=9;        %Blank time (in seconds) will be rounded up to the nearest TR
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

if blanksAfterPass==1;
    totalBlankTime   = length(barOrientations).*ceil(blankTime./params.tr).*params.tr;
    barPassTime=duration.cycle.seconds-totalBlankTime;
    numBarImages=params.numImages-(totalBlankTime./params.tr);
    barPassFrames=barPassTime/duration.stimframe;
else
    barPassTime=duration.cycle.seconds;
    numBarImages=params.numImages;
    barPassFrames=duration.cycle.stimframes;
end



outerRad = params.radius;
innerRad = params.innerRad;
ringWidth = params.ringWidth;

numMotSteps = params.temporal.motionSteps;
numSubRings = params.numSubRings;

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

% r = eccentricity; theta = polar angle
r = sqrt (x.^2  + y.^2);
theta = atan2 (y, x);					% atan2 returns values between -pi and pi
theta(theta<0) = theta(theta<0)+2*pi;	% correct range to be between 0 and 2*pi


% loop over different orientations and make checkerboard
% first define which orientations
orientations = barOrientations./360*(2*pi); % degrees -> rad
remake_xy    = zeros(1,numBarImages)-1;
remake_xy(1:length(remake_xy)/length(orientations):length(remake_xy)) = orientations;
original_x   = x;
original_y   = y;
% step size of the bar
step_nx      = barPassTime./params.tr/length(barOrientations);
step_x       = (2*outerRad) ./ step_nx;
step_startx  = (step_nx-1)./2.*-step_x - (ringWidth./2);
%[0:step_nx-1].*step_x+step_startx+ringWidth./2
disp(sprintf('[%s]:stepsize: %f degrees.',mfilename,step_x));

% Loop that creates the final images
fprintf('[%s]:Creating %d images:',mfilename,numBarImages);
images=zeros(m,n,numBarImages*params.temporal.motionSteps,'uint8');
for imgNum=1:numBarImages
    
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
    
    
    
    loEcc = innerRad;
    hiEcc = outerRad;
    loX   = loX + step_x;
    hiX   = loX + ringWidth;
    
    
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




% Now we compile the actual sequence
% make stimulus sequence, make half and then add the rest as a flipped
% version of the first half
sequence = ...
    ones(barPassFrames./numBarImages,1)*...
    [1:params.temporal.motionSteps:params.temporal.motionSteps*numBarImages];
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
blankCycleFrames=length(sequence)/length(barOrientations);
offTime   = ceil(blankTime./params.tr).*params.tr; % make sure it's a multiple of the tr
offPeriod = ceil(offTime./duration.stimframe);

images(:,:,size(images,3)+1)   = uint8(ones(size(images,1),size(images,2)).*bk);
blankInsert=ones(offPeriod,1).*size(images,3);
newSequence=[];
for passCounter=1:length(barOrientations)
    newSequence=cat(1, newSequence, sequence((passCounter-1)*blankCycleFrames+1:passCounter*blankCycleFrames), blankInsert);
end
sequence=newSequence;

        

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
sequence = [repmat(size(images,3), duration.prescan.stimframes,1); sequence];
timing   = [0:length(sequence)-1]'.*duration.stimframe;
cmap     = params.display.gammaTable;
fixSeq   = [fixSeq(length(fixSeq)+1-duration.prescan.stimframes:end); fixSeq];

% make stimulus structure for output
stimulus = createStimulusStruct(images,cmap,sequence,[],timing,fixSeq);

% save matrix if requested
if ~isempty(params.saveMatrix),
    save(params.saveMatrix,'images');
end;

