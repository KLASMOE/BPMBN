%% testing script for multifactorproductworker() algorithm

% set up factors and indicies:

f1.table = [0 0 0 1 1;...
            0 1 1 0 1;...
            0 0 1 1 1;...
            1 0 0 1 0];
inds1 = [2,4];
v = [0; 1];
f1.logprob = [-1; -2; -4; -8];

f2.table = [0 1 1 1;...
            0 1 0 1;...
            0 1 0 0;...
            1 0 0 0];
inds2 = [1,2];
f2.logprob = [-10; -20; -40; -80];




%% this arranges a new conditional probability table
newtable = [];
% need to construct a multi-column match out of all matching pairs:
shrinktable = f1.table(:,inds1);
matchtable = abs(shrinktable - repmat(v', size(shrinktable,1),1));
rows1 = f1.table(sum(matchtable,2) == 0,:);
p1 = f1.logprob(sum(matchtable,2) == 0);
shrinktable = f2.table(:,inds2);
matchtable = abs(shrinktable - repmat(v', size(shrinktable,1),1));
rows2 = f2.table(sum(matchtable,2) == 0,:);
p2 = f2.logprob(sum(matchtable,2) == 0);
r1 = size(rows1,1);
r2 = size(rows2,1);
rows1 = repmat(rows1,r2,1);
rows2 = repmat(rows2,r1,1);
% complicated indexing achieves staggered repeats for rows1 and rows2
rrconv = repmat([1:r2:r1*r2]',1,r2) + repmat([0:(r2-1)],r1,1);
key = reshape(rrconv,size(rows1,1),1);
rows2 = rows2(key,:);
index1 = true(1,size(rows1,2));
index1(inds1) = false;
index2 = true(1,size(rows2,2));
index2(inds2) = false;
% always put the new column last:
% do this because it is sorted
newtable = [newtable; rows1(:,index1), rows2(:,index2), rows2(:,inds2)];

% multiplication and reshaping actually orders the probabilities in the
% same way as the index convolution above on the tables.
%newp = p1' * p2;
% new logmethod:
if (length(p2) > 1 && length(p1) > 1)
    warning('MultiFactorProductWorker');
end
newp = zeros(length(p1)*length(p2),1);
% multiply (log probs) of each thing in p1 by each thing in p2:
for k = 1:length(p2)
    for j = 1:length(p1)
        newp((k-1)*(length(p1)) + j) = p2(k) + p1(j);
    end
end



%% alternate version:
p1 = p1';
p2 = p2';
newp2 = zeros(size(p1));
srep1 = cell(1,size(p1,2));
for k = 1:size(p1,2)
    repp1 = repmat(p1(:,k),1,size(p2,2));
    srep1{k} = repp1 + p2;
end
for k = 1:size(p1,2)
    for j = 1:size(p2,2)
        newp2(k,j) = logplus(srep1{k}(:,j));
    end
end

colp = reshape(newp2, 1, size(newp2,1) * size(newp2,2));




