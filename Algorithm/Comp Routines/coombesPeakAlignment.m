function [pl, plidx] = coombesPeakAlignment(mz,s, peaks)

%% aligns the peaks in s for easier spike detection...


offset = 20;



N = length(peaks);
pl = zeros(size(peaks));
plidx = pl;



for kk = 1:N
    
    cpeak = peaks(kk);
    
    [~,idx] = min( abs(mz-cpeak)); 
    
    A = s(max(idx-offset-1,1):idx);
    [~,maxi] = max(A);
    
    plidx(kk) = max(idx-offset-2,1)+maxi;
    
    pl(kk) = mz(plidx(kk));
end


% figure(1);
% plot(mz,s); hold on, plot(pl,s(plidx),'o'); plot(peaks,zeros(size(peaks)),'x'),hold off;