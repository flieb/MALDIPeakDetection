function p = detectpeaks_wvlt(s,lambda)
%DETECTPEAKS Peakpicking for MALDI Spectra with WAVELET FRAMES
%   Usage:  [p] = detectpeaks(sig,lambda)
%
%   Input parameters:
%       sig:        Input spectra (rows: spots, columns: spectra for the spot)
%       lambda:     Threshold Parameter
%
%   Output parameters:
%       p:          "Denoised" signal with same size as sig 
%       
%   Description:
%       This method implements the patented algorithm to detect peaks in
%       noisy MALDI data. The output is a "denoised" version of the input
%       signal where the found peaks are marked, the rest of the signal is
%       set to zero. 
%
%   Note:
%       This code has not been vectorized
%
%   See also: 
%
%   Author: F. Lieb, January 2018






%%
[m,n] = size(s);

if (m~=1 && n~=1)
    error('This code only works for a single spectrum');
end


L = length(s);

fmin = 1000;
bw = 1000;
overlap=30;
a = 2;

[g,am,gd,fbr] = waveletframe5('wp2inp',fmin,bw,overlap,a,L,L,1e-6); %was: waveletframe2('wp2inp',fmin,bw,overlap,a,L,L,1e-6,1); 
numwins = (size(g,1)-2)/2;
g = g(2:numwins+1,:);

for kk=1:size(g,1)
    g(kk,:) = g(kk,:)./max(g(kk,:));
end

am1 = ones(size(am));
%plot(abs(g)')

c = cell2mat(wvlttf(s,g',am1)' );

%c = c(:,2:end);

% figure(3);
% ax(1) = subplot(211);
% plot(s);
% ax(2) = subplot(212);
% imagesc(log10(abs(c.'))), set(gca,'YDir','normal');
% linkaxes(ax,'x');


%%
samplesize = 60;
overlap = 0.5;
a = 1;
numintens = L;
M = size(g,1);
numspots = 1;
params.noisest = 0;
params.choice = 0;
%lambda = 3.23e3;

%compute the number of slices needed for signal:
numberofslices = floor( (numintens-samplesize)/floor((1-overlap)*samplesize) + 1 );
idxvec = ones( numberofslices , 2);
idxvec(1,2) = samplesize;
for ii = 2:size(idxvec,1)
    idxvec(ii,:) = idxvec(ii-1,:) + [floor((1-overlap)*samplesize) floor((1-overlap)*samplesize)];
end

%output matrix:
p=zeros(numspots,numintens);
%compute the Gabor transform of the first slice:
lastslf = c(idxvec(1,1) : idxvec(1,2),:);
lastslf = lastslf(:);
nest = 1;

%loop over all slices:
for ii=2:numberofslices %43:43 %2:numberofslices;

    nextslf = c(idxvec(ii,1) : idxvec(ii,2),:);
    nextslf = nextslf(:);
    %edit new
    %nextslice = nextslice - mean(nextslice);
    if params.noisest == 1
        nest = tfnoisest(nextslf,numrows,a,M,samplesize);
    end
    
    m = approxmask(nextslf,lastslf,nest.*lambda,params);

    
    mm = reshape(m,samplesize/a,M);
    %mm = fftshift(abs(mm)-1,1);
    mm = abs(mm)-1;
    c1 = abs(reshape(lastslf,samplesize/a,M))';
    c2 = abs(reshape(nextslf,samplesize/a,M))';
    temp = sum(abs(mm),2);
    c1m = sum(c1);
    c2m = sum(c2);
    differ = c2m - c1m;
    diffsign = sign(differ);
    mul = ones(size(c1m));
    mul(diffsign == -1) = 0;
    peaks = mul'.*temp;
    
    
    p(:,idxvec(ii,1) : idxvec(ii,2)) = peaks';
    
    lastslf = nextslf;
    
   
end



