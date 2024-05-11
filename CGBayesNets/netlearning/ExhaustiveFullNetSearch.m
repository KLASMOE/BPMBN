function [BayesNet, outstats] = ExhaustiveFullNetSearch(contData, discData, ...
    priorPrecision, phencol, searchParameter)
%[BayesNet, outstats] = ExhaustiveFullNetSearch(contData, discData, priorPrecision, phencol, searchParameter)
%
% Performs a greedy hillclimbing search over all possible edges, adding one 
% edges at a time.  In contrast to the usual learning algorithm
% that uses a K2-style node order, this one builds the most likely
% network across all input variables.  Allows backtracking, when removing
% an edge results in an increase in posterior log likelihood of the data.  
% Roughly O(k(n^2)*d*f), where k = priorPrecision.maxParents, n = number
% of nodes, d = size of datasample, and f is the time required for DFS for
% cycle checking, approx ~ O(log(n)).
%
% Input:
% CONTDATA: numeric data matrix of continous nodes (indexed by node,sample)
% DISCDATA: numeric data matrix of discrete nodes (indexed by node,sample)
% PRIORPRECISION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
% PHENCOL: index indicating the phenotype column.  Either the last contdata
%   column or the last discdata column (with contdata before discdata)
% SEARCHPARAMETER: structure with the following fields:
%   searchParameter.BF_THRESH: a Bayes Factor threshold for including an 
%       edge.  Edges with likelihood below BF_THRESH will not be included.  
%       This is the main stopping criterion.
%   searchParameter.nophenotype: a boolean, if true, means that the
%       network should be built without any particular phenotype
%   searchParameter.backtracking: a boolean, if true, means that the
%       network should be built allowing backtracking.  Default = true.
%   searchParameter.DBN: a boolean, if true, indicates that we're learing a
%       two-stage dynamic bayes net, and will use the TSDBNStateTrackerSearch
%       class. Default = false.
%       A dynamic bayes net must further define:
%   searchParameter.t0cont: continuous vars of the first of two-stages
%   searchParameter.tncont: continuous vars of the second of two-stages
%   searchParameter.t0disc: discrete vars of the first of two-stages
%   searchParameter.tndisc: discrete vars of the second of two-stages        
%
% OUTPUT:
% BAYESNET: structure with field:
%   BayesNet.adjMatrix: adjacency matrix for the bayes net.    
% OUTSTATS: structure with fields describing the characteristics of the
%       search procedure, in arrays per "step;" a step is either an edge
%       added or removed.
%   outstats.lldiffs: difference in loglikelihood at each step of the algorithm.
%   outstats.numedges: number of edges in the network at each step of the algorithm.
%   outstats.numevals: number of potential network states evaluated at each
%       step.
%
% Michael McGeachie (c) 2014. MIT license. See cgbayesnets_license.txt.

% phenotype node can't have any edges going in:
if (~isfield(searchParameter,'nophenotype'))
    nophenotype = false;
else
    nophenotype = searchParameter.nophenotype;
end

if (~isfield(searchParameter,'backtracking'))
    backtracking = true;
else
    backtracking = searchParameter.backtracking;
end

if (~isfield(searchParameter,'BF_THRESH'))
    bf_thresh = 0;
    searchParameter.BF_THRESH = 0;
else
    bf_thresh = searchParameter.BF_THRESH;
end

if (~isfield(searchParameter,'DBN'))
    DBN = false;
else
    DBN = searchParameter.DBN;
end

if (~isfield(searchParameter,'unwrapped'))
    unwrapped = false;
else
    unwrapped = searchParameter.unwrapped;
end


% set up the search:
search = StateTrackerSearch(contData, discData, priorPrecision, phencol, backtracking, bf_thresh, nophenotype);
if (DBN)
    if (unwrapped)
        search = UnwrappedTSDBNSTS(contData, discData, priorPrecision, phencol, backtracking, bf_thresh, nophenotype);
        search = search.Init();
        search = search.SetUnwrappedTimeSeries();
    else
        % a Dynamic BN must define:
        % searchParameter.t0cont
        % searchParameter.tncont
        % searchParameter.t0disc
        % searchParameter.tndisc
        % 
        % set up the search to do a dynamic bayesian network:
        search = TSDBNStateTrackerSearch(contData, discData, priorPrecision, phencol, backtracking, bf_thresh, nophenotype);
        search = search.Init();
        search = search.SetTimeSeries(searchParameter.t0cont, searchParameter.t0disc, ...
            searchParameter.tncont, searchParameter.tndisc);
    end
else
    search = search.Init();
end


lldiffs = 0;
numedges = 0;
addedeges = {};
done = false;
count = 0;
evalcount = 0;
evalnums = 0;
%% do an exhaustive search
while (~done && ~search.done)
    count = count + 1;

    % check all possible edges:
    [search, numevals] = search.EvalAllEdges();
    evalcount = evalcount + numevals;
    
    %% decide if an edge should be added:
    lldiffmax = max(max(search.lldiff));
    if (lldiffmax > bf_thresh)
        % adding some edge is a net positive.  Find that edge:
        [~,choseni] = max(max(search.lldiff,[],2));
        [~,chosenj] = max(search.lldiff(choseni,:));
    else
        done = true;
    end

    if (~done)
        [search, success] = search.DoEdgeAdd(choseni, chosenj, lldiffmax);
        if (success)
            lldiffs(end+1) = lldiffs(end) + lldiffmax;
            addedeges{end+1} = [choseni,chosenj];
            numedges(end+1) = search.edges;
            evalnums(end+1) = evalcount;
        end
    end
end


%% set up output to be similar to learnHybridBayesNet()
BayesNet.adjMatrix = search.bn;
BayesNet.weightMatrix = search.weightMatrix;

outstats.lldiffs = lldiffs;
outstats.numedges = numedges;
outstats.numevals = evalnums;
outstats.addededge = addedeges;
outstats.search = search;

