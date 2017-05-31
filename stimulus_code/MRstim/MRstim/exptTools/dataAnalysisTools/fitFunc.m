function [params, q, chisq, df] = fitFunc(params, data, freeParams, fitFuncName, weight)%[params, q, chisq, df] = fitFunc(params, data, freeParams, [fitFuncName], [weight])%% Returns: best-fitting params to the data, a goodness of fit index 'q',% and the chi-square (chisq) and degrees of freedom (df) used to compute q. %% q is the probability that the chisq value of the best-fitting params% is as high as it is due to noise in the data.  In other words, a high% q value (say > .1) means that the model is a good fit.  A low q means% that it is unlikely that the fit error is due to chance, so it's probably% due to a bad model.%% params = starting parameter values%% data = the data to be fit (x in the first row and y in the second)% i.e., data(1,:) = the x values and data(2,:) = the y values.%% freeParams = a list of parameters that are free (the rest are fixed)% (e.g., freeParams = [1 4] will let params(1) and params(4) be free, the% remaining params will be fixed at the value indicated in params).%% fitFuncName = the name of the function to fit.  It defaults to% 'sigmoidLikelihood' if omited or left empty.%% weight = a vector of values used to weight the function fitting.% These can be either a list of standard deviations or the number of% trials at each stimulus level.  Which you send in depends entirely% upon what function fitFuncName is expecting.% (e.g., standard deviations can be used to weight the fit% of certain functions and return a statistically valid chi-squared value% [see Numerical Recipies in C (3rd), pg. 659 - 661])%% An example:% free = [1 2 3];	% list of the free params% fitParams = ones(1, 3)*nan;% slope = 2;% max = 1;% semiSat = mean(x);% startParams = [semiSat slope max];% [fitParams, q, chisq, df] = fitFunc(startParams, [x;y], ...%										free, 'sigmoidLikelihood', sd);%%% For a function to work with this fitting code, it must take the % following form:%%		[Err,chisq,df] = fitFunc(free, data, fixed, [weight])%% where 'free' is a vector of the free parameter values and 'fixed' is% a vector of fixed parameter values with 'nan' used as a placeholder% for the free parameters.  (e.g., the first 'nan' in fixed gets% free(1), the second 'nan' in fixed gets free(2), etc.% data(1,:) is a vector of the x values.% data(2,:) is a vector of the y values%% It should return Err, some measure of the fit error and chisq,% a Chi-square statistic indicating goodness-of-fit, as well% as df, the associated degrees-of-freedom.%% (see sigmoidLikelihood for an example of a function that works with this% fitting code.)%% 99.03.09 RFD: wrote it.% 99.06.15 RFD: moved chi-square calculation to fitFuncName to%					 make this code more general.if ~exist('fitFuncName', 'var')   fitFuncName = [];endif ~exist('weight', 'var')   weight = ones(1,size(data,2));endif isempty(fitFuncName)   fitFuncName = 'sigmoidLikelihood';endif ~isempty(freeParams)   for i=1:length(params)      if any(i==freeParams)         fixed(i) = nan;      else         fixed(i) = params(i);      end   end	params(freeParams) = fminsearch(fitFuncName, params(freeParams), 0, data, fixed, weight);enddf = length(data(1,:)) - length(freeParams);eval(['[L,chisq,df]=' fitFuncName '(params(freeParams), data, fixed, weight);']);% Numerical Recipies in C (3rd ed), pg. 221:% P(observing Chisq by chance for correct model) = 1 - gammainc(chisq/2, df/2)if df < 1   disp('WARNING! fitFunc: you have more freeParams than data points!!!');   df = 1;endq = 1-gammainc(chisq/2, df/2);