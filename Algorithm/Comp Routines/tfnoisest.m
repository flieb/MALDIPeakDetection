function stdn = tfnoisest(c,numrows,a,M,L)

%estimates noise from tf-coefficients c (vector)
%tf coeffs are stored columnwise in c. 

stdn = zeros(1,size(c,2));
for kk = 1:size(c,2)
    A = fftshift(abs(reshape(c(:,kk),L/a,M)),1);
    B = A(1:numrows,:);

    stdn(kk) = median(B(:))/0.6745;
end