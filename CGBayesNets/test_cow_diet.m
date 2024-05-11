% script : Test Cow Diet 'cow_diet.csv'
%
% This script illustrates analysis of the cow_diet sample dataset.  This
% features a non-binary phenotype, which in general I suggest just
% binariziaing; but this will work fine with most of CGBayesNets.
%
% COW DIET is further marked by a large range of values across the
% continuous variables present in the data set.  These need to be
% normalized prior to Bayesian analysis with CGBayesNets.
%

%% MULTIPHENO is set to true to treat the cow diet phenotype as a
% categorical variable with 4 unique values.  set this to FALSE to binarize
% into two categories:
MULTIPHENO = true;

bnpathscript;

[data, cols] = RCSVLoad('cow_diet.csv',false,',');
% first variable, 'Diet' is phenotype
pheno = 'Diet';
% first column is actually a sample ID, drop that:
cols = cols(2:end);

phncol = find(strcmp(pheno, cols));
if (MULTIPHENO)
    % phenotype is coded as : 0, 15, 30, 45.  Needs to be condensed:
    up = unique(data(:,phncol));
    for i = 1:length(up)
        data(data(:,phncol) == up(i)) = i;
    end
else
    % binarize phenotype:
    data(:,phncol) = data(:,phncol) > mean(data(:,phncol));
end

% normalize and center other variables:
for i = 2:length(cols)
    m = mean(data(:,i));
    s = std(data(:,i));
    data(:,i) = (data(:,i) - m) / s;
end

%% then do actual analysis, similar to HybridBN.m:
analysis_title = 'binarized_cow_diet';

fprintf(1,'\n\n======  start new analysis: %s  =====\n\n',analysis_title);
tic;  

% common parameter values:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
priorPrecision.nu = 1;
priorPrecision.sigma2 = 1;
priorPrecision.alpha = 10;
priorPrecision.maxParents = 3;

verbose = true;
doAUC = true;

fprintf(1,'Learning Most Predictive Network Structure for %s\n', analysis_title);
MBNet = LearnStructure(data, cols, pheno, priorPrecision, [analysis_title,'-net']);
fprintf(1,'Learning Network Parameters\n');
MBNet = LearnParamsBN(MBNet);
fprintf(1,'Predicting on Training Data\n');
[acc, p, z] = MBNet.Predict(verbose);

% output AUC or accuracy measure:
% (will output prediction accuracy instead of the AUC for multiple-valued
% phenotypes:
convexhullAUC = AUCWorker(acc,p,z,data(:,phncol),~MULTIPHENO);


