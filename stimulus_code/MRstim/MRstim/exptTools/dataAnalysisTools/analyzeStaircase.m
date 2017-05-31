function analysis = analyzeStaircase(d,varargin);% analysis = analyzeStaircase(d,['doPlot'],['threshErr',numIterations],['fixSlope',slope]);%% Analyzes a single dataSummary from doStaircase (which outputs an array% of dataSummaries).  Fits a Weibull function to the data and returns the% relevant parameters in an analysis struct.  The returned threshold is for% 82% performance.%% Optional arguments can be passed in to effect different actions:%%	'doPlot'	Plots the data and the Weibull fit.  The fit gets returned, too.%%	'threshErr'	Calculates the stderr of the threshold by bootstrapping the data.%				This flag must be followed by the number of iterations with%				which to randomize the data.%% 06.16.99	Written by wap% 07.27.99  RFD: now ignores data points with fewer than 3 trials% dependencies: fitFunc (which calls WeibullLikelihood which calls weib)x = d.stimLevels(:)';k = d.numCorrect(:)';n = d.numTrials(:)';% 07.27.99  RFD: ignore data points with <4 trialsok = find(n>=3); data = [x(ok); k(ok)./n(ok)];% was : thresh = median(x);	% starting value% Changed to median of history to avoid bia to untested levels.% Changed by Junjie Liu May 31, 2002.thresh = median(d.history);	% starting valueguess = 0.5;flake = 0.999;maxThresh = 100;doFixSlope = find(strcmp('fixSlope',varargin));if ~isempty(doFixSlope)	slope = varargin{doFixSlope+1};	free = 1;else	slope = 3.5;	free = [1 2];endstartParams = [thresh slope guess flake maxThresh];[fitParams,q,chisq,df] = fitFunc(startParams, data, free, 'WeibullLikelihood', n(ok));analysis.thresh = fitParams(1);analysis.slope = fitParams(2);analysis.guess = fitParams(3);analysis.flake = fitParams(4);analysis.q = q;analysis.chisq = chisq;analysis.df = df;analysis.x = data(1,:);analysis.y = data(2,:);analysis.n = n(ok);if any(strcmp('doPlot',varargin))	analysis.wx = 10.^[min(log10(data(1,:))):0.01:max(log10(data(1,:)))];	analysis.wy = weib(analysis.wx, analysis.thresh, analysis.slope, guess, flake);	figure;	semilogx(analysis.x,analysis.y,'x');	%lowererr = binomialLowerBound(0.95,round(analysis.y.*analysis.n),analysis.n);	%uppererr = binomialUpperBound(0.95,round(analysis.y.*analysis.n),analysis.n);	%addErrorBars(analysis.x,analysis.y,{analysis.y-lowererr,uppererr-analysis.y});	hold on; plot(analysis.wx,analysis.wy); hold off;enddoThreshErr = find(strcmp('threshErr',varargin));if ~isempty(doThreshErr)	y0 = data(2,:);	n0 = n(ok);	free = 1;	startParams = fitParams;	threshes = [];	fprintf('Doing resampling');	for ii=1:varargin{doThreshErr+1}		for jj=1:length(y0)			data(2,jj) = resample(y0(jj),n0(jj));		end		fitParams = fitFunc(startParams, data, free, 'WeibullLikelihood', n0);		if fitParams(1)<1000			threshes = [threshes fitParams(1)];		end		if ~mod(ii,20)			fprintf('.');		end	end	analysis.threshErr = std(threshes);	fprintf('\n');else	analysis.threshErr = NaN;endreturnfunction y = resample(y0, n)% y0=proportionCorrect, n=numTrials% this function creates an n-long vector with y 1's,% resamples it (with replacement) and returns a new% proportionCorrect valuenCorrect = round(y0*n);% do the easy ones firstif nCorrect==n   y=1;   return;elseif nCorrect==0   y=0;   return;end% create the data vectorz0 = [ones(1,nCorrect), zeros(1,n-nCorrect)];% resample it, with replacementy = 0;for ii=1:length(z0)   y = y+z0(ceil(rand*n));endy = y/n;return