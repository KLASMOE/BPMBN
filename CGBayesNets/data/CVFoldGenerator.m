function [cvdata, cvtest, discvals] = CVFoldGenerator(data, cols, pheno, folds, randfolds, foldsbyclass)
%[cvdata, cvtest, discvals] = CVFoldGenerator(data, cols, pheno, folds, randfolds, foldsbyclass)
%
% Does crossvalidation splitting of data.  
%
% INPUT:
% DATA: data array
% COLS: column names, a cell array of strings
% PHENO: a string for the phenotype, a binary column of DATA, matching one of COLS.
% FOLDS: Number of folds in the cross-validation to perform.  Default = 5.
% RANDFOLDS: boolean, if true, will randomly select each cross validation
%   fold.  not gauranteed to be the same size.  (default = true)
% FOLDSBYCLASS: boolean, if true, will randomize each CV fold within each
% class of the phenotype.
%
% OUTPUT: 
% CVDATA: cell of length FOLDS, contains the training data for fold K in
%   CVDATA{k}.
% CVTEST: cell of length FOLDS, contains the testing data for fold K
% DISCVALS: cell of length COLS, contains an array of values for each
%   discrete column of data
%
%
% Copyright Michael McGeachie, 2016.  MIT license. See cgbayesnets_license.txt.

if (nargin < 4)
    folds = 5;
end
if (nargin < 5)
    randfolds = true;
end
if (nargin < 6)
    foldsbyclass = true;
end


[ncases,~] = size(data);
phcol = find(strcmp(pheno, cols));
phclasses = unique(data(:,phcol));
cases = phclasses(1) == data(:,phcol);
controls = phclasses(2) == data(:,phcol);


mincases = sum(cases) / folds / 2;
mincontrols = sum(controls) / folds / 2;


if (foldsbyclass)
    goodsplit = false;
    while (~goodsplit)
        r1 = ceil(rand(1,sum(cases)) * folds);
        % check that at least minclassmems are in each fold:
        goodsplit = true;
        for i = 1:folds
            sfk = sum(r1 == i);
            if (sfk < mincases)
                goodsplit = false;
                break;
            end
        end
    end
    goodsplit = false;
    while (~goodsplit)
        r2 = ceil(rand(1,sum(controls)) * folds);
        % check that at least minclassmems are in each fold:
        goodsplit = true;
        for i = 1:folds
            sfk = sum(r2 == i);
            if (sfk < mincontrols)
                goodsplit = false;
                break;
            end
        end
    end
    % have two good r1 and r2 vectors; just have to put them together:
    r = -1 * ones(1,ncases);
    r(cases) = r1;
    r(controls) = r2;
else    
    if (randfolds)
        % randomly split data into N-fold CV sets:
        r = ceil(rand(1,ncases) * folds);
    else
        r = zeros(1,ncases);
        for i = 1:folds
            r(ceil(ncases/folds)*(i-1)+1:ceil(ncases/folds)*i) = i;
        end
        r = r(1:ncases);
    end
end

% for each fold in the cross-validation, make sure the discrete nodes have
% the same values:
d = IsDiscrete(data, 5);
discvals = cell(size(d));
for i = 1:length(cols)
    if (d(i))
        % set values of discrete vars:
        discvals{i} = num2cell(unique(data(:,i),'legacy'));
    else 
        discvals{i} = {};
    end
end

cvdata = cell(1,folds);
cvtest = cell(1,folds);
for k = 1:folds
    cvdata{k} = data(r ~= k,:);
    cvtest{k} = data(r == k,:);
end


