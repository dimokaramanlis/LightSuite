function [paramsall, cverr, cvlogli] = fitProbCellModelCv(XX, yy)
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

allpreds = nan(Ndata, 1);
thetaall = nan(Ndata, 1);
for ifold = 1:Ndata
    
    Xtrain = XX(cout.training(ifold), :);
    ytrain = yy(cout.training(ifold), :);
    Xtest  = XX(cout.test(ifold), :);
    % ytest  = yy(cout.test(ifold), :);
    Xtrain = [ones(size(Xtrain, 1), 1) Xtrain];
    Xtest  = [ones(size(Xtest, 1), 1) Xtest];


    pfit = nbreg(Xtrain, ytrain, 'regularization',1e-4);
    allpreds(cout.test(ifold)) = nbvals(pfit, Xtest);
    thetaall(ifold) = pfit.alpha;
    % train
    % pfit = glmfit(Xtrain, ytrain, 'poisson','Options',opts);
    % 
    % % predict
    % allpreds(cout.test(ifold)) = glmval(pfit, Xtest, 'log');
end



paramsall = nbreg([ones(size(XX, 1), 1) XX], yy, 'regularization',1e-4);
alpha    = paramsall.alpha;


mu        = allpreds;  
amu       = alpha*mu;
amup1     = 1 + amu;
cvlogli   = sum( yy.*log( amu./amup1 ) - (1/alpha)*log(amup1) + gammaln(yy + 1/alpha) - gammaln(yy + 1) - gammaln(1/alpha));
cverr     = rsquare(yy, allpreds);
paramsall = paramsall.b;

end