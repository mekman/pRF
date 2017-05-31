function S = dwiComputeSignal(S0, bvecs, bvals, Q)
% Compute expected diffusion signal from tensor
%
%  S = dwiComputeSignal(S0, bvals, bvecs, Q)
%
% This function implements the Stejskal Tanner equation prediction given a
% quadratic form.  There should also be a form of this equation that takes
% in the ADC values, rather than the quadratic form.
%
% This is a version of the Stejskal/Tanner equation for signal attenuation
% See: http://en.wikipedia.org/wiki/Diffusion_MRI#Diffusion_imaging
%
% We need a better description of the expected parameter format (BW).
%
% Parameters
% ----------
% S0:    The signal measured in the non-diffusion weighted scans (B0)  
% bvals: the b values
% bvecs: the b vectors
% Q:     The tensors (quadratic forms) (e.g. see fgTensors) corresponding
% to each node in a voxel. There are often several tensors.  
%
% Returns
% -------
% S: The signal predicted according to this form of the Stejskal/Tanner eq: 
% 
%         S = S0 exp(-bval*(bvec*Q*bvec))
%
%    There is a column of signals for each of the tensors.  So if there are
%    30 directions and 4 tensors, then the returned signals is 30 x 4.
%
% Example:
%
% (c) Stanford VISTA Team, 2012

% Converts the tensors and bvecs into ADC values.  If there are 80
% directions and 4 tensors, the returned ADC is 80 x 4, with each column
% representing the ADCs in all directions for one of the tensors.
% ADC = dtiADC(Q, bvecs);
%
% We have a bval for each ADC:     S = S0 * exp(-bvals .* ADC);
%
% We repmat the bvals to have the same number of rows as Q.  Each row is a
% tensor.  But bvals will have nDirs x nTensors after the repmat.
% S = S0 * exp(- (repmat(bvals, 1, size(Q,1)) .* ADC));
%
S = S0 * exp(- (repmat(bvals, 1, size(Q,1)) .* dtiADC(Q, bvecs)));

% end

end