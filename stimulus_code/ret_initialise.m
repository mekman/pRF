% little script to initialise moving bar stimuli
% PK 18-02-2013

clear all; close all;

% in case there are any bitsis loaded
delete(instrfind)

pRF_paths = fullfile(pwd,'MRstim','MRstim');
pRF_paths = genpath(pRF_paths);
addpath(pRF_paths);

ret