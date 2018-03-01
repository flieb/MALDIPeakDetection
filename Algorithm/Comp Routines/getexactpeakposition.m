function [peaks] = getexactpeakposition(mask,a,M,L,c1,c2)
% takes the approximated mask and extracts the exact position of peaks 
% appearing in the mask, mask consists mainly of zeros, any connected 
% structure apart from zero is considered a peak

%edit 4.11.14 make it work for entire spectrum
%this means mask contains the mask for each spot

[masksize,numberspots] = size(mask);

if (numberspots == 1)

    %since mask is a single row vector, we need to extract 

    mm = reshape(mask,L/a,M);
    mm = fftshift(abs(mm)-1,1);
    c1 = abs(reshape(c1,L/a,M));
    c2 = abs(reshape(c2,L/a,M));

    %figure, subplot(211), imagesc(mm);
    %temp = sum(abs(mm));
    temp = sum(abs(mm));

    c1m = sum(c1);
    c2m = sum(c2);

    %compare sum(c1) with sum(c2) 
    %only consider "new" peaks, meaning it is a peak in c2
    %a big "previous" peak in c1, will cause a false positive in the mask
    %usually with values of opposite sign, to equal out the peak, 
    %we ignore whenever we have such a mask behaviour
    differ = c2m - c1m;
    diffsign = sign(differ);
    mul = ones(size(c1m));
    mul(diffsign == -1) = 0;




    %figure(10), subplot(311), imagesc(c1), subplot(312), imagesc(c2), subplot(313), imagesc(mm);
    %figure(9) , subplot(311), plot(sum(c1)), subplot(312), plot(sum(c2)), subplot(313), plot(mul.*temp);

    %subplot(212), plot(temp);

    %need some more clever idea for edge problems (torus) 
    %causes peaks to be extented to beginning
    %resulting in false peaks... (solved 14.10.14, by ignoring torus like
    %behaviour


    peaks = mul.*temp;
else
    
    mm = reshape(mask,L/a,M,numberspots);
    %mm = fftshift(abs(mm)-1,1);
    %temp = sum(abs(mm));
    temp = sum(abs(fftshift(abs(mm)-1,1)));
    
    %c1 = abs(reshape(c1,L/a,M,numberspots));
    %c2 = abs(reshape(c2,L/a,M,numberspots));
    %c1m = sum(c1);
    %c2m = sum(c2);
    
    %differ = c2m - c1m;
    %diffsign = sign(differ);
    c11 = (reshape(c1,L/a,M,numberspots));
    c22 = (reshape(c2,L/a,M,numberspots));
    %diffsign = sign(sum(abs(c11))-sum(abs(c22)));
    diffsign = sign(sum(abs(c22))-sum(abs(c11)));
    
    mul = ones(size(diffsign));
    mul(diffsign == -1) = 0;
    
    temp2 = mul.*temp;
    
    peaks = squeeze(temp2);


end