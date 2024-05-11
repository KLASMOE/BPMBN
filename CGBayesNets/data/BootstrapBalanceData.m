function balancedData = BootstrapBalanceData(data, cols, pheno)
% balancedData = BootstrapBalanceData(data, cols, pheno)
%
% Take input DATA, which has an unbalanced number of positive/negative
% boolean cases.  It will sample from the underrepresented class with
% replacement to fill out a balanced dataset, with equal numbers of both
% classes.
%
% INPUT:
% DATA: Data matrix, with each sample in rows and each variable in columns
% COLS: cell string of names of each data column
% PHENO: name of the phenotype column, must be a binary column.  Name must
%   match into COLS.
%
% OUTPUT:
% BALANCEDDATA: New datamatrix which is a superset of DATA, with some
% additional rows filled out by sampling from the underrepresented class
%
%
% Copyright Michael McGeachie, 2017.  MIT license. See cgbayesnets_license.txt.


phncol = find(strcmp(pheno,cols));

classes = unique(data(:,phncol));
% not actually required to be binary:
classcount = zeros(size(classes));
for i = 1:length(classes)
    classcount(i) = sum(data(:,phncol) == classes(i));
end

[maxcount, ~] = max(classcount);
newdata = cell(1,length(classes));
for i = 1:length(classes)
    if (classcount(i) < maxcount)
        numnew = maxcount - classcount(i);
%        newdata = zeros(size(numnew,length(cols)));
        cdata = data(classes(i) == data(:,phncol),:);
        r = rand(1,numnew);
        rinds = ceil(r * length(cdata(:,phncol)));
        newdata{i} = cdata(rinds,:);
    end
end

balancedData = data;
for i = 1:length(newdata)
    balancedData = [balancedData;newdata{i}];
end



