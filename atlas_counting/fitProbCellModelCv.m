function [paramsall, cverr] = fitProbCellModelCv(XX, yy)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Ndata      = numel(yy);
% cout       = cvpartition(Ndata, 'KFold', 8);
cout       = cvpartition(Ndata, 'LeaveOut', 1);

accfold    = nan(Ndata, 1);

opts = statset('glmfit');
opts.MaxIter = 250;
warning off;
rng(1);
for ifold = 1:Ndata
    
    Xtrain = XX(cout.training(ifold), :);
    ytrain = yy(cout.training(ifold), :);
    Xtest  = XX(cout.test(ifold), :);
    ytest  = yy(cout.test(ifold), :);


    % train
    pfit = glmfit(Xtrain, ytrain, 'poisson','Options',opts);
    % pfit = glmfit(Xtrain, ytrain, 'binomial','Options',opts);

    % predict
    ypred = glmval(pfit, Xtest, 'log');
    accfold(ifold)    = 2*sum(ypred - ytest + ytest .* log(ytest./ypred));
end
warning on;
paramsall = glmfit(XX, yy, 'poisson','Options',opts);
cverr     = mean(accfold, 'omitmissing');
end