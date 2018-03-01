%% Results Real Data
% this demo script shows how to use the algorithm with spectra-wise peck
% picking as well as how to use the spatially aware approach.


%% 
clear;

% Load Data Set
% must have the following format:
%   - data.x : x-coordinates of spots
%   - data.y : y-coordinates of spots
%   - data.A : NxL Matrix, containing N spectra of length L
%load('ratbrain.mat'); 

lambda = 1.5e-3;

params.noisest=1;


%First, compute the neighors...
coordinates = [data.x data.y, ones(size(ratbrain.x))];
[params.I, D] = rangesearch(coordinates,coordinates,sqrt(2)); %sqrt(2) yields 9 neighbors. 


%Basic Approach: now spatial information included
params.choice = 0;
[p_basic] = detectpeaks_gab(data.A,lambda,60,0.5,15,params);

%Average Filter
params.method = 'average'; %this specifies the method
params.choice = 2;         %this enables spatial awareness
param = 0;
params.wI = getWeights(params.method,D,param);
[p_average] = detectpeaks_gab(data.A,lambda,60,0.5,15,params);

%Gaussian Filter with sigma=0.5
params.method = 'gaussian';
params.wI = getWeights(params.method,D,0.5);
[p_gauss05] = detectpeaks_gab(data.A,lambda,60,0.5,15,params);


%Gaussian Filter with sigma=1
params.method = 'gaussian';
params.wI = getWeights(params.method,D,1);
[p_gauss1] = detectpeaks_gab(data.A,lambda,60,0.5,15,params);

%Median
params.method = 'median';
[p_median] = detectpeaks(data.A,lambda,60,0.5,15,params);


