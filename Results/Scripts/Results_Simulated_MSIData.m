%% Results Simulated MSI Data using Cardinal Package in R
% This m-file reproduces the results for the simulated MSI dataset
% Dataset simulation is based on Cardinal Package in R (see file
% data/simul_MSI.R)


%% create pattern

clear all; clc;
close all

m = 30;
n = 30;


pattern = 0*ones(m,n); %0*

pattern(5:10,4:9) = 2; %2

pattern(19:27,24:24) = 4; %4
pattern(23:23,20:28) = 4; %4

ll = 12;
for kk=1:9
    pattern(14-kk,10+kk:16+ll-kk) = 3; %3
end

lu = 16;
lx = 2;

pattern(lu,10-lx:13-lx) = 1;
pattern(lu+1,8-lx:15-lx) = 1;
pattern(lu+2,7-lx:16-lx) = 1;
pattern(lu+3,6-lx:17-lx) = 1;
pattern(lu+4,6-lx:17-lx) = 1;
pattern(lu+5,5-lx:18-lx) = 1;
pattern(lu+6,5-lx:18-lx) = 1;
pattern(lu+7,5-lx:18-lx) = 1;
pattern(lu+8,6-lx:17-lx) = 1;
pattern(lu+9,6-lx:17-lx) = 1;
pattern(lu+10,7-lx:16-lx) = 1;
pattern(lu+12,10-lx:13-lx) = 1;
pattern(lu+11,8-lx:15-lx) = 1;

clear kk ll lu lx

%% Load Data
%import data
T=readtable('..\..\data\simul_MSI.csv');
data = table2array(T(1:end,2:end));

%peak locations:
vals = [200,400,500,600,800];


%% Visualize data
%visualize ground truth
figure(1)
subplot(441);
imagesc(pattern==1), axis image, axis off;
subplot(442);
imagesc(pattern==2), axis image, axis off;
subplot(443);
imagesc(pattern==3), axis image, axis off;
subplot(444);
imagesc(pattern==4), axis image, axis off;

%visualize noisy data
subplot(445);
imagesc(reshape(data(vals(2),:),m,n)), axis image, axis off;
subplot(446);
imagesc(reshape(data(vals(3),:),m,n)), axis image, axis off;
subplot(447);
imagesc(reshape(data(vals(4),:),m,n)), axis image, axis off;
subplot(448);
imagesc(reshape(data(vals(5),:),m,n)), axis image, axis off;


%visualize some spectra:
figure(2);
subplot(221);
plot(data(:,175)), hold on, plot(vals(2),data(vals(2),175),'*'); hold off;
title('spectrum with peak corresponding to the circle');
subplot(222); 
plot(data(:,156)), hold on, plot(vals(3),data(vals(3),156),'*'); hold off;
title('spectrum with peak corresponding to the square');
subplot(223); 
plot(data(:,520)), hold on, plot(vals(4),data(vals(4),520),'*'); hold off;
title('spectrum with peak corresponding to the triangle');
subplot(224); 
plot(data(:,710)), hold on, plot(vals(5),data(vals(5),710),'*'); hold off;
title('spectrum with peak corresponding to the triangle');

%% do peak picking without spatial awareness

addpath(genpath('..\..\Algorithm'));

%setup parameters:
lambda = 75;
params = struct;
params.choice = 0;
params.noisest = 0;

%perform peak picking
[p] = detectpeaks_gab(data.',lambda,60,0.5,60,params);

%visualize results
figure(1)
subplot(449);
imagesc(reshape(p(:,vals(2)),m,n)), axis image, axis off;
subplot(4,4,10);
imagesc(reshape(p(:,vals(3)),m,n)), axis image, axis off;
subplot(4,4,11);
imagesc(reshape(p(:,vals(4)),m,n)), axis image, axis off;
subplot(4,4,12);
imagesc(reshape(p(:,vals(5)),m,n)), axis image, axis off;

%% do peak picking with spatial awareness

%setup parameters:
lambda = 55;
params = struct;
params.choice = 2;
params.noisest = 0;


%setup neighoring parameters
%get coordinate structure:
coordx = repmat(1:n,m,1);
coordx = coordx(:).';
coordy = repmat(1:m,1,n);
coordz = zeros(1,n*m);
c = [coordx; coordy; coordz].';
[params.I, D] = rangesearch(c,c,sqrt(2)); %sqrt(2) yields 9 neighbors.
params.method = 'gaussian';
param = 0.5;
params.wI = getWeights(params.method,D,param);


%perform peak picking
[p] = detectpeaks_gab(data.',lambda,60,0.5,60,params);

%visualize results
figure(1)
subplot(4,4,13);
imagesc(reshape(p(:,vals(2)),m,n)), axis image, axis off;
subplot(4,4,14);
imagesc(reshape(p(:,vals(3)),m,n)), axis image, axis off;
subplot(4,4,15);
imagesc(reshape(p(:,vals(4)),m,n)), axis image, axis off;
subplot(4,4,16);
imagesc(reshape(p(:,vals(5)),m,n)), axis image, axis off;

