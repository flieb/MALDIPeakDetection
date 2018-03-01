function [G] = stftmat (gin,a,M,L)
%STFTMAT  STFT Analysis Operator with non-torus like behaviour 
%   Usage:  [G] = stftmat(gin,a,M,L)
%
%   Input parameters:
%       gin:    Window function
%       a:      Time sampling step
%       M:      Number of Frequency Samplings
%       L:      Signallength
%
%   Output parameters:
%       G:      Gabor Analysis Operator
%       
%   Description:
%       Produces analysis operator matrix for the Gabor transform. Needed
%       inside detetpeaks(), only usable for small Length sizes, else the
%       operator gets too big.
%       Ignores the torus-effect, since it would result in errors while
%       detecting peaks.
%
%   Example:     
%
%   See also: 
%
%   Author: Florian Lieb, June 2015


%get window from gabwin:
if (iscell(gin))
    g = gabwin(gin,a,M,L);
    windowsize = gin{2};
else
    disp('error: window must be a cell with windowname and size, see gabwin()');
    return;
end

%ensure g is a row vector:
[m,n] = size(g);
if m<n
   g = g.';
end

%remove the part of the signal which is on the right side:
N = max(m,n);
cf = [ones(windowsize,1); zeros(2*N-windowsize,1)];
cf = circshift(cf,-ceil(windowsize/2));
g2 = g.*cf(1:N);

%define the size of output Operator:
l = N/a * M;
b = N/M;
G = zeros(N, l);

%compute G as time- and frequency-shifted versions of g:
for i=1:N/a
    for j=1:M
        G(:, (i-1)*(N/a) + j ) = ifft(circshift(fft(g2),(j-1)*b));
    end
    g = circshift(g,a);
    cf = circshift(cf,a);
    g2 = g.*cf(1:N); %remove the part that gets shifted over the end...
end

end

