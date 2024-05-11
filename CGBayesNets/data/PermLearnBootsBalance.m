function [auc, MBNet, FullBN, outstats] = PermLearnBootsBalance(data, cols, pheno, priorPrecision, verbose, disc, nodelimit)
%[auc, MBNet, FullBN, outstats] = PermLearn(data, cols, pheno, priorPrecision, verbose, disc, nodelimit)
%
% This is an example function for use with the adaptivepermtester() 
% 
% This function performs BN learning after permuting the phenotype
% variable, giving a distribution of possible AUCs for random, scrambled
% data.
%
% INPUT:
% DATA: data array
% COLS: column names, a cell array of strings
% PHENO: a string representing the phenotype column to predict.  Is matched
%   against the COLS array
% PRIORPRECISION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
%       priorPrecision.BFTHRESH: minimum increase in Log Likelihood for 
%           inclusion of an edge in the network.  Deafult = 0;
% VERBOSE: boolean.  If true, increases output.
% DISC: boolean array parallel to COLS, indicating if each column is
%   discrete or not.
% NODELIMIT: hard limit of upper number of nodes allowed to predict the
%   scrambled phenotype
% 
% OUTPUT: 
% AUC: the AUC of this random permutation
% MBNET: the BN, reduced to a markov blanket, identified to predict the 
%   scrambled data 
% MBNET: the full BN identified to predict the scrambled data 
% OUTSTATS: statistics output from the learning algorithm that learned the
%   BN.
%
% Copyright Michael McGeachie, 2016.  MIT license. See cgbayesnets_license.txt.

if (nargin < 5)
    verbose = false;
end
if (nargin < 6)
    nodelimit = 12;
end


% randomize phenotype:
phind = find(strcmp(cols, pheno));

size = nodelimit + 2;
while (size > nodelimit)
    r = randperm(length(data(:,phind)));
    data(:,phind) = data(r,phind);
    bndata = BootstrapBalanceData(data, cols, pheno);
    if (nargin < 6)
        [MBNet, FullBN, outstats] = LearnStructure(bndata, cols, pheno, priorPrecision, '', verbose);
    else
        [MBNet, FullBN, outstats] = LearnStructure(bndata, cols, pheno, priorPrecision, '', verbose, disc);
    end
    size = length(MBNet.mb);
    if (size > nodelimit)
        fprintf('tossing last result, %d is too many nodes!\n', size);
    end
end

MBNet = MBNet.LearnParams();
[acc, p1, z1] = MBNet.Predict(verbose);
auc = AUCWorker(acc, p1, z1, MBNet.GetPhenoCol(), true, false, true, false);
