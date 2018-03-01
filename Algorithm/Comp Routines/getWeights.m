function W = getWeights(filtername,D,param)

%takes filt a 3x3 filter matrix (average, disk, gaussian,...)
%and computes the weights for each neighbor in cell matrix D.

numspots = length(D);


if strcmp(filtername,'average')
    W = arrayfun(@(x) ones(size(D{x}))./numel(D{x}),1:numspots,'UniformOutput',0)';
    return;
end

if strcmp(filtername, 'disk')
    filt = fspecial(filtername,param);
elseif strcmp(filtername, 'gaussian')
    filt = fspecial(filtername,3,param);
else
    error('wrong filter specified');
end

%sort unique filter coeff 
filtunique = sort(unique(filt),'descend');

%compute weight vector according to distance from point to its neighbors
W = arrayfun(@(x) filtunique(ceil(D{x})+1)',1:numspots,'UniformOutput',0)';

            
    