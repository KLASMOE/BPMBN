function tables = MakeDiscContTables(BN)
% tables = MakeDiscContTables(BN)
%
% funciton to loop through the tree of a BN, look at each factor, take the
% root node of that factor, sum out the non-parents of that node, and then
% output the table that remains.
%
% This is intended to compute conditional probability tables for each
% (discrete) node to inspect the posterior likelihood of each variable.
%
% (c) Michael McGeachie 2016.


nodes = BN.nodes;
tree = BN.tree;
tables = cell(size(nodes));
for i = 1:length(nodes)
    % decide if the factor is discrete and what the main node is:
    if (~nodes(i).discrete)
        continue;
    end

    % find all the parents according to the adjmat
    parents = full(BN.adjmat(:,nodes(i).index));
    % also set the self-node to true here so it gets retained:
    parents = logical(parents);
    parents(nodes(i).index) = true;
    pinds = 1:length(BN.adjmat);
    pinds = pinds(parents);

    % else, find a matching factor:
    for j = 1:length(tree)
        missing = setdiff(pinds, tree(j).cpt.findex);
        if (isempty(missing))
            cpt = tree(j).cpt;
        break;
        end
    end

    nonparents = setdiff(cpt.findex, pinds);
    newcpt = cpt;
    for j = 1:length(nonparents)
        % sum out the rest of the factor members using factorsum()
        newcpt = factorsum(newcpt, nodes(nonparents(j)).self);
    end    
    % 
    tables{i} = newcpt;
end


