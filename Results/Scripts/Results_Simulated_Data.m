%% Results Simulated Data from Coombes et al. 
% This m-file creates the results presented in the table, comparing the
% patented method with Wijetunge's algorithm

% The Coombes et al. (2005) data set is available at 
% http://bioinformatics.mdanderson.org/Supplements/Datasets/Simulations/index.html
% in order to recreate the results presented in the dissertation, the first
% 25 data sets have to be used. 


%% Results Wijetunge

clear;

addpath('.\Wijetunge Algorithm');
addpath('.\Comp Routines');


mypath = '.\Data Coombes\';

sensitivity= zeros(2500,1);
fdr = sensitivity;
sigma = 0.5; % as recommended

%Note: Wijetunges algorithm takes 41 hours to process all 2500 spectra!!!
tic
for jj=1:25 % loop over all 25 datasets
    dataset = ['Dataset_' num2str(jj)];
    fprintf('Dataset %2d\n', jj);
    for kk = 1:100 %loop over all 100 spectra
        
        % load spectrum
        numspec = kk; %1-100
        filename = [mypath dataset '\RawSpectra\noisy' num2str(numspec) '.txt'];
        filename2 = [mypath dataset '\truePeaks\truth' num2str(numspec) '.txt'];
        D = table2array(readtable(filename));
        Ds = table2array(readtable(filename2));
        numpeaks = size(Ds,1);
        peaks = Ds(:,1);
        mz = D(:,1);
        s = D(:,2);
        %align peaks to actual maximum:
        [truepeaks,truepeaksidx] = coombesPeakAlignment(mz,s,peaks); 
       
        %Wijetunge's Algorithm:
        pl_wjietjunge3 = PeakDetection(mz,s,sigma);
     
        
        %convert mz values to indices:
        pl_w = zeros(length(pl_wjietjunge3),1);
        for ii=1:length(pl_wjietjunge3)
        [~,pl_w(ii)] = min(abs(pl_wjietjunge3(ii) - mz));
        end

        %compute sensitivity and FDR:
        idx = (jj-1)*100 + kk;
        tmp = evalPeakList(truepeaksidx',40,pl_w');
        sensitivity(idx) = tmp/numpeaks*100;
        fdr(idx) = (length(pl_wjietjunge3) - tmp)./(length(pl_wjietjunge3))*100;

        fprintf('%3d: %5.2f/%5.2f \n',kk,sensitivity(idx),fdr(idx) );
    end
end

%save('wjietunge_SenFdr_0.5_dataset1-25all.mat','sensitivity','fdr');


%% Patent method

clear;

addpath('.\Comp Routines');

mypath = '.\Data Coombes\';


%parameter settings Gabor frame:
Delta = 60;
overlap = 0.5;
wsize = 20;

%general settings
params.choice = 0;  %no neighors
params.noisest = 0; %dont estimate noise


sensitivity= zeros(2500,2);
fdr = sensitivity;

for jj = 1:25
    dataset = ['Dataset_' num2str(jj)];
    fprintf('Dataset %2d\n', jj);
    
    for kk=1:100

        numspec = kk; %1-100
        filename = [mypath dataset '\RawSpectra\noisy' num2str(numspec) '.txt'];
        filename2 = [mypath dataset '\truePeaks\truth' num2str(numspec) '.txt'];
        D = table2array(readtable(filename));
        Ds = table2array(readtable(filename2));
        numpeaks = size(Ds,1);
        peaks = Ds(:,1);
        mz = D(:,1);
        s = D(:,2);
        %align peaks to actual maximum:
        [truepeaks,truepeaksidx] = coombesPeakAlignment(mz,s,peaks); 
        
        %baseline correction (uncomment for baseline corrected data)
        %s = msbackadj(mz,s);

        %detect peaks Gabor frame:
        sig = s./norm(s,1);
        %lambda = 12.8e-9;
        %lambda = 13e-9;
        lambda = 13.1e-9;
        pl_patent = detectpeaks(sig',lambda,Delta,overlap,wsize,params); %wsize was 12
        pl_patent2 = getPeakPositions(pl_patent,sig);
        nlpatent = length(pl_patent2);

        %adjust lambda to number of peaks:
        while (nlpatent<numpeaks)
            if abs(nlpatent-numpeaks)<5
                lambda = lambda*0.85;
            else
                lambda = lambda*0.98;
            end
            [pl_patent] = detectpeaks(sig',lambda,Delta,overlap,wsize,params);
            pl_patent2 = getPeakPositions(pl_patent,sig);
            nlpatent = length(pl_patent2);
            %fprintf('.');
        end

        %detect peaks wavelet frame:
        lambda = 3.1e3;
        %lambda = 3.23e3; %gives slightly better results than 3.5e3
        %lambda = 3.5e3;
        [pl_tmp] = detectpeaks_wvlt(s,lambda);
        pl_pat_wvlt = getPeakPositions(pl_tmp,s);
        nlpatwvlt = length(pl_pat_wvlt);

        while( nlpatwvlt < numpeaks)
            if abs(nlpatwvlt-numpeaks)<5
                lambda = lambda*0.85;
            else
                lambda = lambda*0.98;
            end
            
            [pl_tmp] = detectpeaks_wvlt(s,lambda);
            pl_pat_wvlt = getPeakPositions(pl_tmp,s);
            nlpatwvlt = length(pl_pat_wvlt);
            %fprintf('.');
        end

        %fprintf('Dataset %d: %d Peaks, Gabor: %d, Wvlt: %d\n',numpeaks,length(pl_patent2),length(pl_pat_wvlt));
        
        %compute sensitivity and fdr:
        idx = (jj-1)*100 + kk;
        tmp = evalPeakList(truepeaksidx',40,(pl_patent2)',(pl_pat_wvlt)');
        sensitivity(idx,:) = tmp/numpeaks*100;
        fdr(idx,:) = ([nlpatent nlpatwvlt] - tmp)./([nlpatent nlpatwvlt]).*100;

        %fprintf('Dataset %3d:\n %3d Peaks, Gabor: %3d (S:%2.2f/F:%2.2f), Wvlt: %3d(S:%2.2f/F:%2.2f)\n',kk,numpeaks,length(pl_patent2),sensitivity(kk,1),fdr(kk,1),length(pl_pat_wvlt),sensitivity(kk,2),fdr(kk,2));
        fprintf('%3d: G->%5.2f/%5.2f, W-> %5.2f/%5.2f \n',kk,sensitivity(idx,1),fdr(idx,1),sensitivity(idx,2),fdr(idx,2) );

    end
end



%save('Sensitivity_Fdr_coombesDataset_Patent_PatentWvlt.mat','sensitivity','fdr');
%save('Sensitivity_Fdr_coombesDataset_Patent_PatentWvlt_noBaseline.mat','sensitivity','fdr');



%% computation time of patented algorithm with Gabor frames and vectorized code

clear;

addpath('.\Comp Routines');

%load all 2500 spectra in one matrix:
load('./Data Coombes/ALL_Matrix.mat');

params.choice = 0;
params.noisest = 0;

tic
p = detectpeaks(B,1e-3,60,0.5,20,params);
toc
