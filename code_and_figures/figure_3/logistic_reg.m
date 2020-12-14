% Performs the logistic regression analysis
%--------------------------------------------------------------------------
clear all; close all;

load('all_STP_300ms.mat','Rs','PRRp','sparseness');
nSyn = 20;
all_U = linspace(.05,.9,nSyn);
all_taud = linspace(.02,.5,nSyn);
all_tauf = fliplr(all_taud);

X(:,1) = repelem(all_taud,nSyn*nSyn);
X(:,2) = repmat(repelem(all_tauf,nSyn),1,nSyn);
X(:,3) = repmat(all_U,1,nSyn*nSyn);

SR = [1./sparseness]';
SR(SR>0.5) = 2;
SR(SR<=0.5) = 1;

[B,dev,stats] = mnrfit(X,SR);