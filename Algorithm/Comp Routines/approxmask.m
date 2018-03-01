function [m] = approxmask(c1,c2,lambda,params)
%APPROXMASK  Approximates Gabor Multiplier according to Paper by Doerfler 
%   Usage:  [m] = approxmask(c1,c2,lambda)
%
%   Input parameters:
%       c1:     time-frequency coeffs of current slice
%       c2:     time-frequency coeffs of next slice
%       lambda: thresholding parameter
%
%   Output parameters:
%       m:      Gabor multiplier mask
%       
%   Description:
%       Computes the Gabor multiplier mask according to the paper by
%       Dörfler and Matusiak with
%                   d(m) = || |m| - 1 ||_1
%
%   Example:     
%
%   See also: 
%
%   Author: Florian Lieb, June 2015






if params.choice == 0

% %with full ratbrain this takes 1.2s
% c = abs(conj(c1).*c2);
% %c = abs(temp);
% %angl = angle(temp);
% temp2 = conj(c1).*c1;
% %string = ['1: ' num2str(toc)];
% %disp(string);
% check1 = c./(temp2) >= 1 + ( lambda./(2*temp2) );
% check2 = c./(temp2) <= 1 - ( lambda./(2*temp2) );
% m = ones(size(c1));
% 
% %this takes the most time: 3.4s
% %ignoring the phase speeds it up tremendously, its not really needed:
% val1 = (c - lambda/2)./(temp2);% .*exp(1i*angle(temp));
% val2 = (c + lambda/2)./(temp2);% .*exp(1i*angle(temp));
% 
% %this takes bout 1s 
% m(check1) = val1(check1);
% m(check2) = val2(check2);


%tic
 temp2 = conj(c1).*c1;
 X = abs(conj(c1).*c2)./temp2 -1;
 m = ( X).*max(0,1- (lambda./(2*temp2))./abs( X))+1 ;
%toc

%compnorm(m,m2);

elseif params.choice == 1
    error('wrong choice');
    %structured sparsity in the tf plane, i.e. mask m
    y = (abs(c2.*conj(c1)))./(abs(c1).^2);
    numfn = params.numfn;
    numtn = params.numtn;
    Wc1 = convFreqWeights(c1,numfn,numtn);
    Wc2 = convFreqWeights(c2,numfn,numtn);
    tau2 = lambda./(2*abs(Wc1).^2);
    y2 = Wc2./Wc1;
    m = (y-1).*max(0,1-tau2./(abs(y2-1))) + 1 ; 

elseif params.choice == 2
    %spatial structured sparsity for a single mz value...
    
    %alt1
    %N = arrayfun(@(x) mean(c1(:,params.I(x,:)),2),1:size(c1,2),'UniformOutput',false);
    %M = cell2mat(N);
    
    %alt2
    M = zeros(size(c1));
    for ii=1:size(c1,2)
        %M(:,ii) = mean(c1(:,params.I(ii,:)),2);
        %M(:,ii) = sum(params.wI{ii}.*c1(:,params.I{ii}),2);
        if (strcmp('median', params.method))
            M(:,ii) = median(c1(:,params.I{ii}),2);
        else
            M(:,ii) = sum(params.wI{ii}.*c1(:,params.I{ii}),2);
        end
    end
    
    %alt3
    %M2 = mean(reshape(c1(:,params.I),size(c1,1),size(c1,2),size(params.I,2)),3);
    %compnorm(M2,M);
    
    X = abs(conj(c1).*c2)./(conj(c1).*c1);
    tX = conj(M).*M;
    cX = abs(conj(M).*c2);
    
    m = (X -1).*max(0,1-(lambda./(2*tX.*abs(cX./tX-1))))+1 ;

    
    % ----
%      temp2 = conj(c1).*c1;
%      X = abs(conj(c1).*c2)./temp2 -1;
%      m2 = ( X).*max(0,1- (lambda./(2*temp2))./abs( X))+1 ;
%      compnorm(m,m2);
    
end

end