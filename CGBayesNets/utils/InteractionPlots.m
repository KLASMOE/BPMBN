function topnvals = InteractionPlots(topnvals, snpdata, rsnames, phdata, phenotitle, repsnps, reppheno)
% INPUT:
% 
% TOPNVALS: structure with fields:
%   topnvals.vals = interaction strength values
%   topnvals.rows
%   topnvals.cols

if (nargin > 5)
    topnvals.reprows = repsnps(:,topnvals.rows);
    topnvals.repcols = repsnps(:,topnvals.cols);
end


totalsnps = size(snpdata,2);
totalp = totalsnps^2 /2;

% filter down to just the good ones:
topnvals.datarows = snpdata(:,topnvals.rows);
topnvals.datacols = snpdata(:,topnvals.cols);
topnvals.rownames = rsnames(topnvals.rows);
topnvals.colnames = rsnames(topnvals.cols);
topnvals.pheno = phdata;


priorPrecision.nu = 10;
priorPrecision.sigma2 = 1;
priorPrecision.alpha = 10;
priorPrecision.maxParents = 3;

% check the bayesian score for these:
[int_best, main1, main2, reverse_best] = BayesianInteraction(topnvals.pheno, ...
    topnvals.datarows, topnvals.datacols, priorPrecision);

% if there is replication SNP data, check the bayesian score from THOSE:
if (nargin > 5)
    int_rep = BayesianInteraction(reppheno, topnvals.reprows, ...
        topnvals.repcols, priorPrecision);
    topnvals.int_rep = int_rep;
end

topnvals.int_best = int_best;
topnvals.reverse_best = reverse_best;

% do a histogram and normal-tail plot:
figure();
k1 = length(topnvals.vals) / totalp;
x1 = min(topnvals.vals);
k2 = length(topnvals.vals(1:100))/ totalp;
x2 = min(topnvals.vals(1:100));
q1 = norminv(1-k1,0,1);
q2 = norminv(1-k2,0,1);
sigma = (x2 - x1)/(q2 - q1);
mu = x1 - q1 *((x2-x1)/(q2-q1));
h = histogram(topnvals.vals);
hold on;
dist = max(topnvals.vals) - min(topnvals.vals);
y = [min(topnvals.vals) : dist/1000 : max(topnvals.vals)];
f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(y,f * h.BinWidth * totalp,'LineWidth',1.5);
hold off;
legend('Interaction','Normal Fit');
title([phenotitle, ' Top Interactions (n = ',num2str(length(topnvals.vals)),' of ',num2str(totalp),')']);
xlabel('Interaction Strength (MI)');
ylabel('count');

% convert to p-vals using that distribution:
pvals = 1-normcdf(topnvals.vals,mu,sigma);
pvalqqplot(pvals, true, [], false, totalp);
title([phenotitle, ' Interaction Strength p-values QQ-plot']);

% plot all interactions by main effects and by BDE
figure();
plot(topnvals.vals,int_best,'ro');
xlabel('Interaction Strength');
ylabel('Bayesian Dirichlet Equivalent (Uniform) log Likelihood Diff');
title([phenotitle, ' Interaction Effects']);


%figure();
%plot(topnvals.vals,reverse_best,'bo');
%xlabel('Interaction Strength');
%ylabel('Bayesian Dirichlet Equivalent (Uniform) log Likelihood Diff');
%title('Reverse Interaction Effects');

figure();
% add min here so this is >= 0.
valscaling = 15 + 100 * (topnvals.vals.^2 - min(topnvals.vals.^2)) / (max(topnvals.vals.^2) - min(topnvals.vals.^2));
scatter(main1, main2, valscaling, int_best);
colorbar;
xlabel('Main Effect 1 (BDEu)');
ylabel('Main Effect 2 (BDEu)');
title([phenotitle, 'Bayesian Interaction (Color, BDEu); MI Interaction (Size, MI)']);

if (nargin > 5)
    figure();
    plot(topnvals.int_best,topnvals.int_rep,'ro');
    xlabel('Discovery Bayesian Interaction Strength');
    ylabel('Replication Bayesian Interaction Strength');
    title([phenotitle, ' Bayesian Replication']);
end


% regresion test:
topnvals.b = zeros(size(topnvals.rownames));
topnvals.brep = zeros(size(topnvals.rownames));
for i = 1:length(topnvals.rownames)
    y = topnvals.pheno;
    
    X = [ones(size(topnvals.datarows(:,i))),topnvals.datarows(:,i),topnvals.datacols(:,i),...
        topnvals.datarows(:,i) .* topnvals.datacols(:,i)];
    [b,bint] = regress(double(y),double(X));
    topnvals.b(i) = b(end); % last one is the interaction effect
    
    if (nargin > 5)
        y = reppheno;

        X = [ones(size(topnvals.reprows(:,i))),topnvals.reprows(:,i),topnvals.repcols(:,i),...
            topnvals.reprows(:,i) .* topnvals.repcols(:,i)];
        [b,bint] = regress(double(y),double(X));
        topnvals.brep(i) = b(end); % last one is the interaction effect
    end
    
end

% plot regression vals vs. MI vals:
figure();
plot(topnvals.vals,topnvals.b,'ro');
xlabel('Interaction Strength (MI)');
ylabel('Regression Test Interaction (beta)');
title([phenotitle, ' Interaction Effects']);

if (nargin > 5)
    figure();
    plot(topnvals.b,topnvals.brep,'ro');
    xlabel('Discovery Interaction (Regression beta)');
    ylabel('Replication Interaction (Regression beta)');
    title([phenotitle, ' Replication of Interaction Effects']);    
end


if (nargin > 5)
    figure();
    plot(topnvals.vals,topnvals.int_rep,'ro');
    xlabel('Discovery Interaction (MI) ');
    ylabel('Replication Interaction (BDE-u)');
    title([phenotitle, ' Replication of Interaction Effects']);    
end




