 function params = setRetinotopyParams(expName, params)
% setRetinotopyParams - set parameters for different retinotopy scans 
%
% params = setRetinotopyParams([expName], [params])
%
% Sets parameter values for the specified expName.  
%
% params is a struct with at least the following fields:
%  period, numCycles, tr, interleaves, framePeriod, startScan, prescanDuration
%
% Returns the parameter values in the struct params.
% If called with no arguments, params will be a cell array listing
% all the experiment names that it is configured to do.
%
% 99.08.12 RFD rewrote WAP's code with a cleaner wrapper.
% 05.07.04 SOD ported to OSX; several changes
% 10.02.19 BMH Modified for Ress lab, cut out unnecessay parts

% the following should match those listed in the switch statement below

expNames = {'rotating wedge (45deg duty)','rotating wedge (90deg duty)','rotating wedge (45deg duty) Reverse', ...
            'rotating wedge with blanks (45deg duty)','rotating wedge with blanks (45deg duty) Reverse', ...
			'expanding ring (45% duty)','expanding ring (90% duty)', 'contracting ring (45% duty)', ...
			'expanding ring with blanks (45% duty)', 'contracting ring with blanks (45% duty)', ...
			'full-field, on-off', 'full-field, on-off (impulse)',...
            'Translating Bars 1', 'Translating Bars 2', 'Translating Bars 3', 'Translating Bars 8 Pass (Dumoulin)'};

if ~exist('expName', 'var')
	params = expNames;
	return;
end
disp(['[' mfilename ']:Setting stimulus parameters for ' expName '.']);

% some more definitions
if isfinite(params.interleaves),
    params.framePeriod = params.tr*params.interleaves;
else,
    params.framePeriod = params.tr;
end;
params.startScan   = 0;
params.quitProgKey = KbName('q');


%disp('flipping images to simulate 3T projector view');
%params.flipUD = 1;
%params.flipLR = 1;

%--------------------------------------
% background id, you can change this for manual calibration when only
% having three intensities (black, white and gray=bg)
bg = 128;
%bg = 144;
%--------------------------------------

if ~isempty(params.calibration),
    params.display = loadDisplayParams('displayName',params.calibration);
    disp(sprintf('[%s]:loading calibration from: %s.',mfilename,params.calibration));
else,
    params.display.screenNumber   = max(Screen('screens'));
    [width, height]=Screen('WindowSize',params.display.screenNumber);
    params.display.numPixels  = [width height];
    params.display.dimensions = [24.6 18.3];
    params.display.pixelSize  = min(params.display.dimensions./params.display.numPixels);
    params.display.distance   = 43.0474;%40;
    params.display.frameRate  = 60;
    params.display.cmapDepth  =  8;
    params.display.gammaTable = [0:255]'./255*[1 1 1];
    params.display.gamma      = params.display.gammaTable;
    params.display.backColorRgb   = [bg bg bg 255];
    params.display.textColorRgb   = [255 255 255 255];
    params.display.backColorRgb   = bg;
    params.display.backColorIndex = bg;
    params.display.maxRgbValue    = 255;
    params.display.stimRgbRange   = [0 255];
    params.display.bitsPerPixel   = 32;
    disp(sprintf('[%s]:no calibration.',mfilename));    
end;
params.display.quitProgKey = params.quitProgKey;

if max(Screen('screens')) < params.display.screenNumber,
    disp(sprintf('[%s]:resetting screenNumber %d -> %d.',mfilename,...
        params.display.screenNumber,max(Screen('screens'))));
    params.display.screenNumber   = max(Screen('screens'));
end;

% IMPORTANT: Set stereoFlag to 1 if using stereo display.  This     %
% will affect both the stimulus presentation and the fixation point %
params.stereoFlag = 0;
params.display.stereoFlag = 0;


% Flickering fixation point parameters
%
% this controls the duration of the fix flicks, in frames.
% Set it to 0 to turn get no flicks.
params.fixFlickFrames = 5;
% this controls the density of the flicks, per frame.
% Thus, .01 will, on average, flick once every 100 frames.
params.fixFlickFreq = .01;

params.dispString = [expName '.  Please watch the fixation square.'];

%
% Color parameters
%
params.backRGB.dir = [1 1 1]';	% These two values are your
params.backRGB.scale = 0.5;		% standard default gray.
params.stimLMS.dir = [1 1 1]';
params.stimLMS.scale = 1.0;
%bk = findName(params.display.reservedColor,'background');
%params.display.reservedColor(bk).gunVal = (params.display.numColors-1) * ...
%								params.backRGB.scale*params.backRGB.dir';

%
% Set some common defaults

params.temporal.frequency = 2; %Hz
params.temporal.motionSteps = 10;
if ischar(params.stimSize),
    if isfield(params.display, 'Rect')
        params.radius = pix2angle(params.display, floor(min((params.display.Rect(3)-params.display.Rect(1))/2, (params.display.Rect(4)-params.display.Rect(2))/2))); 
    else
        params.radius = pix2angle(params.display,floor(min(params.display.numPixels)/2));
    end
else
    params.radius = params.stimSize;	
end;

if isfield(params.display, 'Rect')
    tmpy=(params.display.Rect(4)+params.display.Rect(2))/2;
    tmpx=(params.display.Rect(3)+params.display.Rect(1))/2;
    tmpsize=angle2pix(params.display, params.radius);
    params.display.Rect=[tmpx-tmpsize tmpy-tmpsize tmpx+tmpsize tmpy+tmpsize];
    params.display.Rect=params.display.Rect(:);
end

disp(sprintf('[%s]: Stimulus size: %.1f degrees / %d pixels.',...
    mfilename,params.radius,angle2pix(params.display,2*params.radius)));
% front projector=16; back projection=3;
params.seqDirection = 0;	% 0 or 1- just plays the movie backwards if set to 1

% Wedge parameters
params.innerRad = 0;		% Non-zero for annular wedge condition (deg)
params.wedgeDeg = 90;		% Wedge polar angle (deg)
params.subWedgeDeg = 15;	% Sub wedge polar angle (deg) 

% Ring parameter - 8 for a radius=16 stim gives a 180 degree duty cycle
params.ringDeg = params.radius/2;			% Ring radius/width (deg)

% Wedge and ring parameters
params.subRingDeg = 1;			% 1/2 radial spatial freq (deg)

params.dynamicStaticFlag = 0;	% if set, half the time will be a static stim

params.numImages = params.period/params.framePeriod;  % Number of samples of the image (i.e. per cycle)
params.duration = params.period/params.numImages;

% insert blanks
params.insertBlanks.do = 0;
params.insertBlanks.freq = 4;
params.insertBlanks.phaseLock = 0;

switch expName
case 'rotating wedge (90deg duty)',
	params.type = 'wedge';		% Set to 'wedge' or 'ring'
	params.wedgeDeg = 90;
	params.seqDirection = 0;
    params.numSubRings = (params.radius)/(2*params.subRingDeg);
case 'rotating wedge (90deg duty) Reverse',
	params.type = 'wedge';		% Set to 'wedge' or 'ring'
	params.wedgeDeg = 90;
	params.seqDirection = 1;
    params.numSubRings = (params.radius)/(2*params.subRingDeg);
case 'rotating wedge (45deg duty)',
	params.type = 'wedge';
	params.wedgeDeg = 45;
	params.seqDirection = 0;
	params.innerRad = 0;
    params.numSubRings = (params.radius)/(2*params.subRingDeg);
case 'rotating wedge (45deg duty) Reverse',
	params.type = 'wedge';
	params.wedgeDeg = 45;
	params.seqDirection = 1;
	params.innerRad = 0;
    params.numSubRings = (params.radius)/(2*params.subRingDeg);
case 'rotating wedge with blanks (45deg duty)',
	params.type = 'wedge';
	params.wedgeDeg = 45;
	params.seqDirection = 0;
	params.innerRad = 0;
    params.insertBlanks.do = 1;
    params.numSubRings = (params.radius)/(2*params.subRingDeg);
case 'rotating wedge with blanks (45deg duty) Reverse',
	params.type = 'wedge';
	params.wedgeDeg = 45;
	params.seqDirection = 1;
	params.innerRad = 0;
    params.insertBlanks.do = 1;
    params.numSubRings = (params.radius)/(2*params.subRingDeg);
case 'expanding ring (90% duty)',
	params.type = 'ring';
	params.ringDeg = params.radius/4;
	params.seqDirection = 0;
    params.numSubRings = (params.radius)/(2*params.subRingDeg);
case 'contracting ring (180deg duty)',
	params.type = 'ring';
	params.ringDeg = params.radius/2;
	params.seqDirection = 1;
    params.numSubRings = (params.radius)/(2*params.subRingDeg);
case 'expanding ring (45% duty)',
	params.type = 'ring';
	params.ringDeg = params.radius/8;
	params.seqDirection = 0;
    params.numSubRings = (params.radius-params.innerRad)/(params.radius);
case 'contracting ring (45% duty)',
	params.type = 'ring';
	params.ringDeg = params.radius/8;
	params.seqDirection = 1;
    params.numSubRings = (params.radius-params.innerRad)/(params.radius);
case 'expanding ring with blanks (45% duty)',
	params.type = 'ring';
	params.ringDeg = params.radius/8;
	params.seqDirection = 0;
    params.insertBlanks.do = 1;
    params.numSubRings = (params.radius-params.innerRad)/(params.radius);
case 'contracting ring with blanks (45% duty)',
	params.type = 'ring';
	params.ringDeg = params.radius/8;
	params.seqDirection = 1;
    params.insertBlanks.do = 1;
    params.numSubRings = (params.radius-params.innerRad)/(params.radius);
case 'full-field, on-off',
	params.type = 'center-surround';
	params.centerInnerRad = 0;
	params.centerOuterRad = params.radius;
	params.surroundInnerRad = params.radius;
	params.surroundOuterRad = params.radius;
	params.numImages = 2;
	params.duration = params.period/params.numImages;
    params.numSubRings = (params.radius)/(2*params.subRingDeg);
case 'full-field, on-off (impulse)',
	params.type = 'center-surround';
	params.centerInnerRad = 0;
	params.centerOuterRad = params.radius;
	params.surroundInnerRad = params.radius;
	params.surroundOuterRad = params.radius;
	params.numImages = 2;
	params.duration  = params.period/params.numImages;
    params.numSubRings = (params.radius)/(2*params.subRingDeg);
    params.impulse   = input('Please enter impulse time (sec): ');
    params.insertBlanks.do = 1;    
    params.insertBlanks.freq = params.numCycles;
    params.insertBlanks.phaseLock = 1;
case {'Translating Bars 1', 'Translating Bars 2', 'Translating Bars 3'},
	params.type = 'bar';
	params.ringDeg = 2; 
    params.seqDirection = 0;
    params.insertBlanks.do = 1;
    params.numSubRings = (params.radius-params.innerRad)/(params.radius);
case {'Translating Bars 8 Pass (Dumoulin)'},
    params.type = 'bar';
	params.ringDeg = params.radius./4; 
	params.seqDirection = 0;
    params.insertBlanks.do = 1;
    params.numSubRings = (params.radius-params.innerRad)/(params.radius);   
otherwise,
	error('Unknown expName!');
end


% stimulus on/off presentation
if params.insertBlanks.do,
%    bn = questdlg('Phase lock stimulus on/off cycle?','Blanking','Yes','No','No');
%	if strmatch(bn,'Yes'),
%		params.insertBlanks.phaseLock = 1;
%    else,
		params.insertBlanks.phaseLock = 0;
%    end;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%																	%
% Fixation parameters												%
%																	%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.display.fixType        = params.fixation; % disk or largeCross
params.display.fixSizePixels  = 6;%3;%6;12
switch(lower(params.display.fixType))
    case {'disk','double disk'}
        params.display.fixColorRgb    = [255 0 0 255;...
                                         0 255 0 255];%172 0 0  255];

        dim.x = params.display.numPixels(1);
        dim.y = params.display.numPixels(2);
        params.display.fixX = round(dim.x./2);
        params.display.fixY = round(dim.y./2);
        
       if isfield(params.display, 'Rect')
            params.display.fixY=round((params.display.Rect(2)+params.display.Rect(4))/2);
        end

    case {'small cross +'}
        params.display.fixColorRgb    = [255 0 0 255;...
                                         0 255 0 255];
        dim.x = params.display.numPixels(1);
        dim.y = params.display.numPixels(2);
        params.display.fixX = round(dim.x./2);
        params.display.fixY = round(dim.y./2);
        if isfield(params.display, 'Rect')
            params.display.fixY=round((params.display.Rect(2)+params.display.Rect(4))/2);
        end       
        

    otherwise,
        error('Unknown fixationType!');
end



params.fix.task               = 'Detect fixation change';
params.fix.colorRgb           = params.display.fixColorRgb;
params.fix.responseTime       = [.01 3]; % seconds


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%																	%
% Calculations (not to be updated by user)							%
%																	%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.ringWidth=params.ringDeg;

% Polar angle of wedge in radians
params.wedgeWidth = params.wedgeDeg * (pi/180);

% Number of rings in each wedge
%params.numSubRings = (params.radius-params.innerRad)/(2*params.subRingDeg);
%params.numSubRings = (params.radius-params.innerRad)/(params.radius);
%params.numSubRings = (params.radius)/(2*params.subRingDeg);

% Number of wedges in each ring
params.numSubWedges = params.wedgeDeg/(2*params.subWedgeDeg);

% duration of each image (seconds) 
params.imageDuration = params.period / params.numImages; 

% Duration of params
params.scanDuration = params.period * params.numCycles + params.prescanDuration;

% some checks, must be done before we reset certain params
retParamsCheck(params);

% ***HACK!  We'll let makeRetinotopy add the prescan stuff
params.ncycles = params.numCycles;
params.prescanDuration = params.prescanDuration;
params.period = params.period;


