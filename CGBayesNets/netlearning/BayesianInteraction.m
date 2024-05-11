function [int_best, main1, main2, reverse_best] = BayesianInteraction(pheno, ...
    childData1, childData2, priorPrecision)
% [int_best, main1, main2, reverse_best] = BayesianInteraction(pheno, ...
%    childData1, childData2, priorPrecision)
%
% Compute interaction (Bayesian Posterior Likelihood) for two children over
% and above the main effects on the phenotype.
%
% INPUT:
%   PHENO: column of phenotype data (discrete)
%   CHILDDATA1: array of possible child variables
%   CHILDDATA2: array of other pairs of child variables.  Each pair of
%   (CHILDDATA1(:,i), CHILDDATA2(:,i)) will be checked for interaction
%   PRIORPRECISION: optional bayesian parameters.
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%
%
% OUTPUT:
%   INT_BEST: the interaction value above the main effects of the pair of
%       variables.  INT_BEST(i) is the interaction value of (CHILDDATA1(:,i), 
%       CHILDDATA2(:,i)).
%
% Copyright Michael McGeachie, 2013.  MIT license. See cgbayesnets_license.txt.




if (nargin < 4)
    priorPrecision.nu = 10;
    priorPrecision.sigma2 = 1;
    priorPrecision.alpha = 10;
end

main1 = zeros(size(childData1(1,:)));
main2 = main1;
int_best = main1;
reverse_main1 = main1;
reverse_main2 = main1;
reverse_best = main1;

for i = 1:length(childData1(1,:))
    % compute baseline likelihood
    cData1 = double(childData1(:,i));
    cData2 = double(childData2(:,i));
    model_ini_child1=learnLocalBN_DiscToDisc([],cData1',priorPrecision);
    model_ini_child2=learnLocalBN_DiscToDisc([],cData2',priorPrecision);


    % compute main effects
    parentData = double(pheno);
    model_addedge1=learnLocalBN_DiscToDisc(parentData',cData1',priorPrecision);
    model_addedge2=learnLocalBN_DiscToDisc(parentData',cData2',priorPrecision);

    % main effects:
    main1(i) = model_addedge1.logLLH - model_ini_child1.logLLH;
    main2(i) = model_addedge2.logLLH - model_ini_child2.logLLH;


    % compute interaction BDE:
    model_int1=learnLocalBN_DiscToDisc([parentData,cData2]',cData1',priorPrecision);
    model_int2=learnLocalBN_DiscToDisc([parentData,cData1]',cData2',priorPrecision);

    int1 = model_int1.logLLH - model_ini_child1.logLLH;
    int2 = model_int2.logLLH - model_ini_child2.logLLH;

    int_best(i) = max(int1,int2);

    % do reverse (phenotype is child)
    model_ini_reverse=learnLocalBN_DiscToDisc([],parentData',priorPrecision);

    % reverse main effects:
    model_reverse1=learnLocalBN_DiscToDisc(cData1',parentData',priorPrecision);
    model_reverse2=learnLocalBN_DiscToDisc(cData2',parentData',priorPrecision);
    reverse_main1(i) = model_reverse1.logLLH - model_ini_reverse.logLLH;
    reverse_main2(i) = model_reverse2.logLLH - model_ini_reverse.logLLH;

    % reverse interactions:
    reverse_int=learnLocalBN_DiscToDisc([cData1,cData2]',parentData',priorPrecision);
    reverse_int1 = reverse_int.logLLH - model_reverse1.logLLH;
    reverse_int2 = reverse_int.logLLH - model_reverse2.logLLH;
    reverse_best(i) = max(reverse_int1,reverse_int2);
end
