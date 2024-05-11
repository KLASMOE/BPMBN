function tables = GetProbTables(BN, fid)
%PrintCPTs(BN, fid)
%
% Take a hybrid bayesian network and get the conditional probability tables
% for all the nodes.  Output them in a simple cell array.
%
% INPUT:
% BN: A BayesNet class object representing the bayes net.  Must be
%   completed with LearnParamsBN().
% FID: (optional) a file identifier from FOPEN(FNAME), to output the tables
%   to.  Otherwise, will print tables to standard out.
%
% OUTPUT:
% TABLES: the conditional probability tables; or list of LPPotentials for 
%   continuous nodes.
%
% (c) Michael McGeachie 2016.

if (nargin < 2)
    fid = 1;
end

disctables = MakeDiscContTables(BN);
nodes = BN.nodes;
tree = BN.tree;

% have to loop through the nodes, and find each node somewhere in the tree
tables = cell(size(nodes));
for i = 1:length(nodes)
    if (nodes(i).discrete)
        tables{i} = disctables{i};
        printdisc(tables{i}, nodes(i), fid);
    else
        for j = 1:length(tree)
            % ClusterSetTree objects
            if ((tree(j).index == i) || ...
                (tree(j).discrete && sum (tree(j).dmembers == i) > 0))
                % found the node
                % get the cpt:
                tables{i} = tree(j).lppotential;
                printcont(tables{i},nodes,fid);
                break;
            end
        end
    end
end
end

function printdisc(table, node, fid)
% assume that the table.table is empty
% print a headerline:
probs = exp(table.logprob);
probs = probs / sum(probs);
fprintf(fid,'Table for Node %s:\n',node.self);
for i = 1:length(table.factors)
    fprintf(fid,'%s\t',table.factors{i});
end
fprintf(fid,'p()\n');
% loop through all assignments:
[v, vinds] = ValueIncrement(table.values, []);
while (~isempty(v))
    % convert the combinatorial values assignment into an index:
    factorinds = 1:length(v);
    spind = multisparsefactorinds(table, factorinds, v);
    for j = 1:length(table.factors)
        fprintf(fid,'%d\t',v(j));
    end
    p = full(probs(spind));
    fprintf(fid,'%f\n',p);

    [v, vinds] = ValueIncrement(table.values, v, vinds);
end

end


function printcont(table, nodes, fid)
for i = 1:length(table)
    % table is a list of lppotentials, for each value of conditioning
    % discrete vars:
    lpp = table{i};
    fprintf(fid,'p(%s',nodes(lpp.head).self);
    if (~isempty(lpp.conditionvars))
        fprintf(fid,' | ');
        for j = 1:length(lpp.conditionvars)
            fprintf(fid, '%s = %d', nodes(lpp.conditionvars(j)).self, ...
                nodes(lpp.conditionvars(j)).values{lpp.conditionvalinds});
            if (j < length(lpp.conditionvars))
                fprintf(fid, ', ');
            end
        end
    end
    fprintf(fid,') ~ N(');
    fprintf(fid, '%f', lpp.const);
    for j = 1:length(lpp.tail)
        fprintf(fid, ' + %f %s', lpp.params(j), nodes(lpp.tail(j)).self);
    end
    fprintf(fid, ', %f)\n', lpp.sigma);
end
end




