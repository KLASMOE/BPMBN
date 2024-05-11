% test script for missing data (by blocks)
% work for Rose Du

%% IMPORTANT: 
%   check if we have to normalize LLH computations by the number of
%   elements, which would have always been constant under an assumption of
%   complete data
%   ANSWER: probably not?

% CHECK : anywhere function unique() is used, each NaN adds a separate
% element


%% load MissingData:
% try it as a txt
% [numdata, cols] = RCSVLoad('MissingData.txt',false,'\t');
% this isn't going to work with leading missingness, since the textscan()
% function trims leading whitespace

% so far, have fixed (first pass) :
% 1. RCSVLoad() : now accepts missing values in CSV input
% 2. IsDiscrete() : now skips NaN's so they don't auto-set columns to
%       continuous
% 3. learnLocalBN_ContToCont() : Actually just sets NaNs to 0, which then
%       drop out of the matrix multiplies.  This solution saves information
%       on the rest of a row with some NaNs in it.  (the other option being
%       to drop any row with >=1 NaN value.)
% 4. learnLocalBN_DiscToDisc() : No Change; calls LLMultiParent()
% 5. LLMultiParent() : Changed so NaN values will match with anything;
%       making them count as all possible values in rows.  This preseverse
%       information in a row that contains a NaN value (the other option
%       being to just drop a row if it has >=1 NaN value).
%       ALSO changed so that it returns 0 (no information) if >= 1 parent
%       has 100% missing values, which should prevent adding extranous
%       parents with high missingness.
% 5. initDiscreteFactor() : didn't seem to need editing; worked just right
%       by ignoring NaNs anyway
% 6. PushEvidence() : changed to just exit if the evidence is empty (NaN).
% 7. PredictPheno() : changed to use the new function MultiFactorSum() to
%       sum out all the remaining discrete vars that had no evidence
%       entered on them. 
% 8. InitContFactor() : changed to not compute means of data without first
%       filtering out NaNs.
% 9. InitDiscreteFactor() : Seemed to be an error that assumed there would
%       be only 1 variable in the factor, which normalized successive pairs
%       of logprobs in the CPT.  Changed to a global normalization across
%       the whole CPT, all logprobs.
% 10. LPPotential.MakeWeight() : changed to throw an error if the input
%       value is NaN.
% 11. LPPotential.AddEvidence() : changed to throw an error if the input
%       value is NaN.
% 12. HybridDiscreteEvidence() : changed to just do nothing if all the
%       input is NaN or empty.  Will propagate partially empty input.
% 13. factorreduceworker() : now drops out all NaN evidence before
%       proceding.
% 14. ClusterSetTree() : unchanged
% 15. MultiFactorSum() : NEW FUNCTION.  Repeatedly calls factorsum(), which
%       isn't new, but wasn't used.
% 16. multifactorproductworker() : Changed to deal with larger clusters of
%       factors, which seems to be a consequence of missing data.  Updated 
%       core algorithms. 
% 17. getDiscParentConfig() : Changed "unique()" invocation to not
%       consider NaNs, which add spurious elements into the result of the 
%       unique() call 
% 


%% try it as a csv:
[numdata, cols] = RCSVLoad('MissingData_v3.csv',false,',',false,true);


pheno = cols{1};
data = numdata;
alg = 1;
BFTHRESH = 0;
verbose = true;

[auc,MBNet,BNet] = BNLearn(data, cols, pheno, alg, BFTHRESH, verbose);

%% Check normal Binary Network for comparison:
[numdata, cols] = RCSVLoad('BinaryVars.csv',false,',',false,true);

pheno = cols{1};
data = numdata;
alg = 1;
BFTHRESH = 0;
verbose = true;

[auc,MBNet,BNet] = BNLearn(data, cols, pheno, alg, BFTHRESH, verbose);



%% Check large files from Rose Du

fname = '../../privatedata/rose/aneu_test_v1.csv';
[numdata, cols, strdata, numcolinds, strcolinds] = RCSVLoad(fname,false,',',false,true);

pheno = cols{1};
data = numdata;
alg = 1;
BFTHRESH = 0;
verbose = true;

[auc,MBNet,BNet] = BNLearn(data, cols, pheno, alg, BFTHRESH, verbose);



%% Check larger file from Rose Du
% this one is 6000 features

fname = '../../privatedata/rose/aneu_au_blood_pheno_rnaseqp05_rrbsp05_intersect_sep_tidnan.csv';
[numdata, cols, strdata, numcolinds, strcolinds] = RCSVLoad(fname,false,',',false,true);

pheno = cols{1};
data = numdata;
alg = 1;
BFTHRESH = 0;
verbose = true;

[auc,MBNet,BNet] = BNLearn(data, cols, pheno, alg, BFTHRESH, verbose);


%%
% this one is 650 features

fname = '../../privatedata/rose/aneu_au_blood_rnaseqp05_rrbsp05_intersect_septtid_z.csv';
[numdata, cols, strdata, numcolinds, strcolinds] = RCSVLoad(fname,false,',',false,true);

pheno = cols{1};
data = numdata;
alg = 1;
BFTHRESH = 0;
verbose = true;

[auc,MBNet,BNet] = BNLearn(data, cols, pheno, alg, BFTHRESH, verbose);


