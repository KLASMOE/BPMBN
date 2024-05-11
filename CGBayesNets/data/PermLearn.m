function [auc, MBNet, FullBN, outstats] = PermLearn(data, cols, pheno, priorPrecision, verbose, disc, nodelimit, algorithm, filter)
%[auc, MBNet, FullBN, outstats] = PermLearn(data, cols, pheno, priorPrecision, verbose, disc, nodelimit, algorithm, filter)
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
    disc = IsDiscrete(data);
end
if (nargin < 7)
    nodelimit = 12;
end
if (nargin < 8)
    algorithm = 1;
end
if (nargin < 9)
    filter = false;
end

outstats = [];

% randomize phenotype:
phind = find(strcmp(cols, pheno));
phcol = data(:,phind);

size = nodelimit + 2;
while (size > nodelimit)
    r = randperm(length(phcol));
    data(:,phind) = data(r,phind);
    
    % if filtering, filter here :
    if (filter)
        bfs = BayesFactorScore(data, cols, pheno, priorPrecision, disc);
        bfs(phind) = inf;
        [~,bfsort] = sort(bfs,'descend');
        maxnum = min(length(cols),2*nodelimit);
        bndata = data(:,bfsort(1:maxnum));
        bncols = cols(bfsort(1:maxnum));        
    else
        bndata = data;
        bncols = cols;
    end
    
    % use one of 5 network learning algorithms:
    if (algorithm == 1)
        [MBNet, FullBN, outstats] = LearnStructure(bndata, bncols, pheno, priorPrecision, '', verbose, disc);
    end
    if (algorithm == 2)
        [MBNet, FullBN, outstats] = LearnPhenoCentric(bndata, bncols, pheno, priorPrecision, priorPrecision.BFTHRESH, verbose, disc);
    end
    if (algorithm == 3)
        [FullBN, outstats] = FullBNLearn(bndata, bncols, pheno, priorPrecision.BFTHRESH, '', priorPrecision, disc, verbose);
        MBNet = FullBN.MakeIntoMB();
    end
    if (algorithm == 4)
        MBNet = NBBayesNet(bndata, bncols, pheno, priorPrecision.BFTHRESH, '', priorPrecision, disc);    
        FullBN = MBNet;
    end
    if (algorithm == 5)
        MBNet = TreeAugmentedNB(bndata, bncols, pheno, priorPrecision.BFTHRESH, '', priorPrecision, disc);
        FullBN = MBNet;
    end
    size = length(MBNet.mb);
    if (size > nodelimit)
        fprintf('tossing last result, %d is too many nodes!\n', size);
    end
end




MBNet = MBNet.LearnParams();
[acc, p1, z1] = MBNet.Predict(verbose);
auc = AUCWorker(acc, p1, z1, MBNet.GetPhenoCol(), true, false, true, false);
