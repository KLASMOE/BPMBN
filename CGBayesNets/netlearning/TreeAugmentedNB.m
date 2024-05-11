function TreeNBNet = TreeAugmentedNB(data, cols, pheno, BFTHRESH, outfilename, priorPrecision, disc, verbose, varlim)
%TreeNBNet = TreeAugmentedNB(data, cols, pheno, BFTHRESH, outfilename, priorPrecision, disc, verbose, varlim)
%
% 
% function to build a Tree Augmented Naive Bayes model from data
%
% WARNING: Doesn't allow edges between Disc and Cont nodes
% Needs to be fixed to work on mixed networks, by not suggesting edges
% that are from continuous nodes to discrete nodes.
%
% Copyright Michael McGeachie, 2020.  MIT license. See cgbayesnets_license.txt.
% 

if (nargin < 7)
    disc = IsDiscrete(data);
end
if (nargin < 8)
    verbose = false;
end
if (nargin < 9)
    varlim = inf;
end


% first get a NB model to start from:
NBNet = NBBayesNet(data, cols, pheno, BFTHRESH, outfilename, ...
    priorPrecision, disc, verbose, varlim);


phncol = find(strcmp(pheno,NBNet.cols));

% then take all the elements that aren't the phenotype,
% compute LLH for added an additional parent:
llh = -inf * ones(size(NBNet.adjmat));
for i = 1:length(NBNet.cols)
    for j = 1:length(NBNet.cols)
        % test out i as the parent of j
        if (i == j)
            continue;
        end

        % already using the phenotype column:
        if (i == phncol || j == phncol)
            continue;
        end
        % can't have a continuous parent and discrete child
        if (NBNet.disc(j) && ~NBNet.disc(i))
            continue;
        end
        % And for now: also prevent disc->cont edges
        if (NBNet.disc(i) && ~NBNet.disc(j))
            continue;
        end
        
        
        % and the phenotype is already the parent of j
        childData = NBNet.data(:,j);       
        contParentData = [];
        discParentData = NBNet.data(:,phncol);
        % get baseline LLH:
        % compute base LLH
        if (~NBNet.disc(j))
            % j is a continuous node
            model_base=learnLocalBN_MixToCont(contParentData',discParentData',childData',NBNet.priorPrecision);
        else
            model_base=learnLocalBN_DiscToDisc(discParentData',childData',NBNet.priorPrecision);
        end
        
        
        if (NBNet.disc(i))
            discParentData = [discParentData,NBNet.data(:,i)];
        else
            contParentData = NBNet.data(:,i);
        end

        % compute new LLH with the proposed edge:
        if (~NBNet.disc(j))
            % j is a continuous node
            model_new=learnLocalBN_MixToCont(contParentData',discParentData',childData',NBNet.priorPrecision);
        else
            model_new=learnLocalBN_DiscToDisc(discParentData',childData',NBNet.priorPrecision);
        end
        % save the bayes factor here, which is the improvement in
        % likelihood:
        llh(i,j) = model_new.logLLH - model_base.logLLH;
    end
end

llhorig = llh;

% drop phenotype:
phncol = strcmp(pheno,NBNet.cols);
llh = llh(~phncol,~phncol);
% then, we negate this
llh = -1 * llh;
% and maybe convert to an undirected graph by averaging values?
for i = 1:length(llh)
    for j = i+1:length(llh)
        if (~isinf(llh(i,j)) && ~isinf(llh(j,i)))
            avval = mean([llh(i,j),llh(j,i)]);
            llh(i,j) = avval;
            llh(j,i) = avval;
        end
    end
end

% need to fix this later, probably build our own min span tree alg, which
% works on directed graphs.

gllh = graph(llh,NBNet.cols(~phncol));
% and use the MATLAB minimum spanning tree algorithm:
% this function requires a symmetric matrix
[mintree, parents] = minspantree(gllh);

% then take the topology identified here and convert it into a tree ontop
% of what was already known in the NBNet:
TreeNBNet = NBNet;

% need to walk through the PARENTS and convert to directed edges agian on
% the NBN:
phnind = find(phncol);
for i = 1:length(parents)
    % add one to indices that appear after the phenotype node was removed
    if (parents(i) == 0)
        continue;
    end
    if (phnind > i)
        nodeind = i;
    else
        nodeind = i + 1;
    end
    if (phnind > parents(i))
        parentind = parents(i);
    else
        parentind = parents(i) + 1;
    end
    % now add the edge:
    TreeNBNet.adjmat(parentind,nodeind) = 1;
    TreeNBNet.weightMatrix(parentind,nodeind) = llhorig(parentind, nodeind);
    
end

% need to rebuild the Tree and the nodes arrays in the BN:
TreeNBNet = TreeNBNet.RebuildNodes();
TreeNBNet.tree = [];  % just make sure the tree is dead, will get rebuilt JIT





