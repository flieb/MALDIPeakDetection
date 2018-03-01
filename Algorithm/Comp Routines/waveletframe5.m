function [g,am,gd,fbr,lfb,ufb] = waveletframe5(wname,fmin,bw,overlap,a,L,fs,epsilon)

%generates the waveletframe and the dual frame

if nargin<8
    epsilon = 1e-6;
end

x = linspace(1/fs,L/fs,L);
fmin = fmin*L/fs^2;
bw = bw*L/fs^2;

unscaled = 0;

[wfun,k1,kend,j] = mwavelet(wname,a,fmin,bw,L/fs,1e-6);

step = (1/(j*overlap));

tmp = k1: step : kend;
scales = length(tmp);
g = zeros(2*scales+2,length(x));
am = ones(2*scales + 2,1);
idx = 2;
for k = tmp
    y = wfun(x,k,j);
    y = real(y);
    y(isnan(y)) = 0;
    [~,~,s] = support(y,epsilon);
    am(idx) = helperfun(s,L);
    am(end-idx+2) = am(idx);
    g(idx,:) = sqrt(am(idx)).*y;
    %g(idx,:) = y;
    g(end-idx+2,:) = fliplr(g(idx,:));
    idx = idx + 1;
end

temp = bsxfun(@times,g,sqrt(1./am));

p = sum(temp.^2);

[~,r] = support(g(2,:));
indx1 = findnextmax(r+1,p);
am(1) = helperfun(2*indx1,L);
l = support(g(scales+1,:));
indx2 = findnextmax(l-1,p);
am(scales+2) = helperfun(L-2*indx2,L);

%temp = bsxfun(@times,g,sqrt(1./am));

%p = sum(temp.^2);
tmp = max(p) - p;

g(1,1:indx1) = sqrt(am(1).*tmp(1:indx1));%tmp(1:indx1);%
g(1,end-indx1+1:end) = sqrt(am(1).*tmp(end-indx1+1:end));%tmp(end-indx1+1:end);%sqrt(am(1).*tmp(end-indx1+1:end));
g(scales+2,indx2:length(x)-indx2+1) = sqrt(am(scales+2).*tmp(indx2:length(x)-indx2+1));%tmp(indx2:length(x)-indx2+1);%sqrt(am(scales+2).*tmp(indx2:length(x)-indx2+1));

% p = sum(g.^2);
% 
% for ii=2:scales+1
%     am(ii) = estimateAM(g(ii,:),epsilon * min(p)/(2*scales + 2));
%     am(end-ii+2) = am(ii);
% end
%     
% g = bsxfun(@times,g,sqrt(am));

temp = bsxfun(@times,g,sqrt(1./am));

p = sum(temp.^2);
%figure(4), plot(p); pause()

fbr = max(p)-min(p);
lfb = min(p);
ufb = max(p);
gd = bsxfun(@times,g,1./p);

if unscaled
    
    g = bsxfun(@times,g,1./sqrt(am));
end



function y = helperfun(x,L)
%finds the correct a<s such that a is divisible by L
xl = [1 alldiv(L)];
indx = find( xl-L/x <= 0,1,'last');
if isempty(indx)
    y = 1;
else
    y = xl(indx);
end


function indx = findnextmax(idxin,p)
%finds the next maximum such that the two cover functions are smoother

indx = idxin;
%looking for the next max
if (p(indx) <= p(indx+1))
    %rising flank
    while(p(indx) < p(indx+1)), indx = indx + 1; end
else
    %falling flank
    while(p(indx) < p(indx-1)), indx = indx - 1; end
end

    

function am = estimateAM(g,epsilon)
%finds the correct am such that ||I-S||_2 is leq eps
[r,c] = size(g);
[~,posmax] = max(g);
if (r>c), g = g.'; end
if (norm(g) ~= 1), g = g/norm(g); end
g = abs(g).^2;
tmp2 = cumsum((g));
il = find(tmp2 >= epsilon,1,'first');

if isempty(il)
    am = length(g);
else
    deltax = posmax - il;
    am = helperfun(deltax,length(g));
end

