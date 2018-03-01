function [p,nest] = detectpeaks_gab(sig,lambda,samplesize,overlap,wsize,params)
%DETECTPEAKS Peakpicking for MALDI Spectra with GABOR FRAMES
%   Usage:  [p] = detectpeaks(sig,lambda)
%           [p] = detectpeaks(sig,lambda,samplesize,overlap,wsize)
%
%   Input parameters:
%       sig:        Input spectra (rows: spots, columns: spectra for the spot)
%       lambda:     Threshold Parameter
%
%    (optional:)
%       samplesize: Length of slices
%       overlap:    Overlap of slices
%       wsize:      Size of window in Gabor Transform
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
%   Example:     
%       
%   See also: 
%
%   Author: F. Lieb, January 2018


if nargin < 3
    samplesize = 60;
    overlap = 0.5;
    wsize = 20;
    params.choice = 0;
    params.noisest = 0;
end

if nargin == 5
    params.choice = 0;
    params.noisest = 0;
end


%get number of spots and number of spectra
[numspots, numintens] = size(sig);
if (numintens == 1)
    sig = sig.';
end


%Settings for the Gabor slices to compute the Multiplier from
a = 1;
M =  samplesize;
g = {'hann',wsize};

%compute Analysis Matrix for Gabor Frame G(g,a,M):
FM = stftmat(g,a,M,samplesize); 

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
firstslice = sig(:,idxvec(1,1) : idxvec(1,2));
lastslf = FM'*firstslice';

if (numspots ~= 1)
  h = waitbar(0,'Please wait', 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
  setappdata(h,'canceling',0)
end




%noise estimation from tf-cofficients
numrows = round(0.1*samplesize);
nest = 1;

%loop over all slices:
for ii=2:numberofslices %43:43 %2:numberofslices;

    nextslice = sig(:,idxvec(ii,1) : idxvec(ii,2));
    nextslf = FM'*nextslice' + eps;
    if params.noisest == 1
        nest = tfnoisest(nextslf,numrows,a,M,samplesize);
        if nest<100*eps
            nest = 1;
        end
    end
    
    m = approxmask(nextslf,lastslf,nest.*lambda,params);
 
    temp2 = getexactpeakposition(m,a,M,samplesize, lastslf,nextslf);
    p(:,idxvec(ii,1) : idxvec(ii,2)) = temp2';
    
    lastslf = nextslf;
    
    if (numspots ~= 1)
        if getappdata(h,'canceling')
             break
        end
        string = ['slice ' num2str(ii) '/' num2str(numberofslices)];
        waitbar(ii/numberofslices,h, string);
    end
end
if (numspots ~= 1)
    delete(h);
end

end


