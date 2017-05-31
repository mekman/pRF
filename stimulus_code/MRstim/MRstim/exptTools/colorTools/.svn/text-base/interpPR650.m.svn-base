function [measure, wavelengths] = interpPR650(measureRaw)
% 
% [measure, wavelengths] = interpPR650(measureRaw)
% 
%AUTHOR: Poirson
%DATE:   11.07.97
% 
% This function converts the raw PR-650 measurements into the form we
% conventionally use at Stanford for our color calculations.  The
% interpolated values are in watts/sr/m^2/nm.  
% 
%   To compute candelas/m2 from the returned vector, measure, use the formula
% 
%     L (cd/m^2) = 683 * Integral [ V(lambda)*Radiance(lambda) dLambda]
% 
%  or, in matlab terms, candm2 = 683 * vlambda' * measure
% 
% Why is there a 683?  Because 1 watt/sr/m^2 at the peak of
% vLambda (555nm) is defined to be 683 cd/m^2 by the CIE.
% 
%ARGUMENTS:
%   INPUT:
%   measureRaw:   Data produced by the PR-650.  This is a Nx2 matrix.  The first
%                 column contains the wavelength positions and the second column
%                 contains the measured values in watts/sr/m^2.  The wavelengths
%                 always run from 380-780 in 4nm steps.
% 
%   RETURNS:
%   
%   measure:     These are interpolated data
%   wavelengths: These are the wavelengths, so that
%      plot(wavelengths,measure) is the result
%   
% 
%CHANGES:
% 
% 09.09.98 (BW):  I updated the values to maintain the correct
% physical units.  The returns were made column vectors.
% 
%   Discussion:
% 
% The values the PR 650 file are measured in watts/sr/m^2/4nm
% band  The relationship between this quantity and candelas/m^2
% is 
% 
%     L (cd/m^2) = 683 * Integral [ V(lambda)*Radiance(lambda) dLambda]
% 
%  with radiance measured in watts/sr/m^2/lambda
% 
% Since the data are in 4nm wide bands, when we interpolate to 1
% nm bands we should divide by 4.  Then we multiply by 683.  The
% combination of these two operations is approximately 170 ~
% 683/4. 
% 
% 12.01.97 (ABP) You might want to multiply measurements by 170
% so that dot product with XYZ produces proper photometric values.
% I don't know why this seems to be the magic scale factor.
% 

disp('InterpPR650 was updated to use the correct units on Sept. 9, 1998');
disp('See Wandell if you have any questions.')
  
% PR-650 data start in range [380-780] in increments of 4nm.  The
% units are watts/sr/m^2 per 4nm band.  Our goal is to create a
% file that runs from 370 to 730 for consistency with our
% software.  This might be a mistake, but there you have it.

prWave	= measureRaw(:,1);	% Measured wavelengths in nm
prValues= measureRaw(:,2);	% Measured values in watts/sr/m^2

% linearly interpolate over the range that intersects the 
% PR range (380:780) and what we use (370:730) to get 1nm spacing
% 
w1nm 	= [380:730];
mInterp	= interp1(prWave,prValues,w1nm)/4;

% Zero pad front of interpolated spectrum
% 
zeroPadFront = zeros(10,1)';
waveFront = [370:379];

wavelengths = [waveFront,w1nm]';
measure     = [zeroPadFront, mInterp]';

return
