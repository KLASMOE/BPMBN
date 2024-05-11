function [localModel]=learnLocalBN_ContToCont(parentData,childData,priorPrecision)
%[localModel]=learnLocalBN_ContToCont(parentData,childData,priorPrecision)
%
% This function estimates parameters and computes log-likelihood of the 
% regression-like Bayesian network, given the data. The network assumes a 
% structure where multiple/no Gaussian parent(s) modulate a single 
% Gaussian child.
%
% INTERNAL
%
% INPUT:
% PARENTDATA: numeric data matrix of data for the parents of the node
% CHILDDATA: numeric data matrix of data for the one child node
% PRIORPRECISION: a structure including the usual HybridBayesNets
%   parameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % hard-limit on the number of parents
%           each node
%
% OUTPUT:
%   localModel.logLLH: log likelihood of the model with PARENTS connected 
%       to CHILD;
%   localModel.postMean: [];
%   localModel.postVar: posterior variance of the model
%   localModel.postPrec: tau; 
%   localModel.regreCoeff: beta_new;
%%
% Copyright Hsun-Hsien Chang, 2010.  MIT license. See cgbayesnets_license.txt.



%% check input data
numChild = size(childData,1);
numSample = size(childData,2);
numParent = size(parentData,1);
numSampleP = size(parentData,2);

if (numSampleP>0) && (numSample~=numSampleP)
    error('The numbers of samples in parent and child data are not same.');
else
    clear numSampleP;
end

if numChild ~= 1
    error('There must be 1 child node.');
end





%% check priorPrecision
useSampleStat = false; %true=sample statitics; false=Bayesian method
if isempty(priorPrecision) 
    useSampleStat = true;
else 
    if isempty(priorPrecision.nu) && isempty(priorPrecision.sigma2)
        useSampleStat = true;
        nu0=1;
        sigma2_0=1;
    elseif isempty(priorPrecision.nu) && ~isempty(priorPrecision.sigma2)
        nu0=1;
        sigma2_0 = priorPrecision.sigma2;
    elseif ~isempty(priorPrecision.nu) && isempty(priorPrecision.sigma2)
        nu0=priorPrecision.nu;
        sigma2_0=1;
    else
        nu0=priorPrecision.nu;
        sigma2_0=priorPrecision.sigma2;
    end
end




%% arrange the data
% first drop rows of parent data that have NaNs:
X = [ ones(numSample,1), parentData'];
y = childData';

% Just replace NaN with 0's, they drop out of the equation anyway
X(isnan(X)) = 0;
y(isnan(y)) = 0;


TESTING = false;
if (useSampleStat || TESTING)
    %% regression coefficient
    beta_new = regress(y,X); % actually this skips NaN's just fine

    %% conditional mean 
    tau = mean(X,1)*beta_new;
    
    %% conditional variance
    yEst=X*beta_new; % estimates of y
    sigma2_new =  sum((y-tau).^2)/(numSample-1);
            
    %% compute log-likelihood (which is log-probability)
    % aka: -0.5*numSample (1 + log(2pi*sigma2_new))
    logLLH = -0.5*log(2*pi*sigma2_new)*numSample - 0.5*sum((y-tau).^2)/sigma2_new;
    logLLHdataest = logLLH;
end

    
if (~useSampleStat) % use Bayesian approach
    % equataions from Ferrazzi, Sebastiani, Ramoni, Bellazzi, "Bayesian 
    % approaches to reverse engineer cellular systems: a simulation study 
    % on nonlinear Gaussian networks." BMC Bioinformatics, 2007.
    

    %% prior parameters
    beta0 = zeros(numParent+1,1); %% regression coefficients 
    R0 = eye(numParent+1);  
    alfa1 = nu0/2;
    alfa2 = 2/(nu0*sigma2_0); 

    %% update parameters
    alfa1n = nu0/2 + numSample/2;
    Rn = R0 + X'*X;
    % previously was inv(Rn) * (R0*beta0+X'*y), but MATLAB claims this is
    % better:
    betan = Rn\(R0*beta0+X'*y);
    inv_alfa2n = (-betan'*Rn*betan+y'*y+beta0'*R0*beta0)/2 + 1/alfa2;
    alfa2n = 1/inv_alfa2n;
    nun = nu0 + numSample;
    sigma2_n = 2/(nun*alfa2n);

    %% Bayesian estimates of parameters
    tau = alfa1n*alfa2n;
    beta_new = betan;
    sigma2_new = nun*sigma2_n/(nun-2);

    %% compute log-likelihood
    %Prob=1/((2*pi)^(n/2))*(sqrt(det(R0))/sqrt(det(Rn)))*(gamma(nun/2)/gamma(nu0/2))*((nu0*sigma0/2)^(nu0/2))/((nun*sigman/2)^(nun/2));
    logLLH = -numSample/2*log(2*pi) + ...
              0.5*( log(det(R0))-log(det(Rn)) ) + ...
              gammaln(nun/2) - gammaln(nu0/2) + ...
              nu0/2*log(nu0*sigma2_0/2) - nun/2*log(nun*sigma2_n/2);
%    fprintf('LLH data: %f, LLH Bayes: %f\n', logLLHdataest, logLLH);
end


%% output
localModel.logLLH = logLLH;
localModel.postMean = [];
localModel.postVar = sigma2_new;
localModel.postPrec = tau; 
localModel.regreCoeff = beta_new;
