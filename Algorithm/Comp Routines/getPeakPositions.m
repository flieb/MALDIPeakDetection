function pl = getPeakPositions(in,s)

%in: Input from patent
%s: orig signal
%pl: peak list



tmp = in;
tmp(in > 0) = 1;


idx = find(tmp);

if (sum(idx) == sum(1:length(s)) )
    pl = 0;
    warning('peak position not possible, peaks not seperated by zeros');
    return ;
end

idx = idx(diff([1 idx]) ~= 1);

pl = zeros(length(idx),1);

for kk = 1:length(idx)
    
    %find endindx
    begidx = idx(kk);
    endidx = 0;
    while (tmp(begidx + endidx) == 1 && begidx+endidx < length(s))
        endidx = endidx + 1;
    end
    
    %get maximum from s
    A = s(begidx:begidx+endidx);
    [~,midx] = max(A);
    
    pl(kk) = begidx+midx-1;
end