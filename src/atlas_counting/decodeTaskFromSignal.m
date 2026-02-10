function [decacc] = decodeTaskFromSignal(Xkin, ychoice)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%----------------------------
%%
[Nareas, Nmice] = size(Xkin);
Ncps = 4;
%-------------------------------------------------------------------------
Xdata       = zscore(Xkin');
[~, pcvec]  = pca(Xdata);
%-------------------------------------------------------------------------
rng(1);
ikeep  = abs(ychoice) > 0;
cout   = cvpartition(nnz(ikeep), "Leaveout");
% cout   = cvpartition(nnz(ikeep), 'KFold', 5);

lambdas = logspace(-6,2,10);

xinput = pcvec(ikeep, 1:Ncps);
mdl = fitclinear(xinput, ychoice(ikeep, :) >0,...
    'Learner', 'logistic', 'CVPartition',cout, 'Verbose',0,...
    'ObservationsIn','rows', 'Solver','sparsa', 'Lambda', lambdas);

[decacc, imax] = max(1-kfoldLoss(mdl));




end