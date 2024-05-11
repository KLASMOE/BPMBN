function [gif, cout, sorder] = pvalqqplot(pvals, logs, colorscale, nogif, totaln, fdrthresh, pointlabels, forcepointlabels)
% [gif] = pvalqqplot(pvals, logs, colorscale, nogif, totaln)
% 
% do a log-log QQ plot of p-values vs. a uniform distribution
%
% INPUT:
% PVALS: list of pvalues
% LOGS: boolean, should be "true", plots negative log pvalues.  Can leave
%   blank.
% COLORSCALE: will plot x's based on color specified by this value;
%   parallel to PVALS
% NOGIF: boolean, by defaut FALSE.  If true, will not print the GIF and
%   LAMBDA on the graph.
% TOTALN: (optional) the input PVALS represent the top length(PVALS) of the
%   total number of tests run.  This plots just the upper tail of the
%   qqplot.
% FDRTHRESH: (optional) pvalues representing tests passing an FDR (False
%   Discovry Rate) threshold of less than FDRTHERSH will be highlighted with
%   a red circle. (DEFAULT = 0.05).
%
% OUTPUT: 
% GIF: Genome-Inflation Factor.  This is a measure of the stratification
%   within the population; or the general systematic enrichment of the
%   p-values for significance.  If (much) greater than 1.0, it means the null
%   distribution does not hold; and therefor the modeling assumptions are
%   incorrect.
%
% Copyright Michael McGeachie, 2012.  MIT license. See cgbayesnets_license.txt.


if (nargin < 2 || isempty(logs))
    logs = true;
end
if (nargin < 3 || isempty(colorscale))
    docolor = false;
else
    docolor = true;
end
if (nargin < 4 || isempty(nogif))
    nogif = false;
end
if (nargin < 5 || isempty(totaln))
    totaln = length(pvals);
end
if (nargin < 6 || isempty(fdrthresh))
    fdrthresh = 0.05;
end
if (nargin < 7)
    pointlabels = {};
end
if (nargin < 8)
    forcepointlabels = {};
end


% flip input
if (size(pvals,1) ~= 1)
    pvals = pvals';
end

if (logs)
    sample = (1:length(pvals)) ./ (totaln +1);
    x = -1*log10(sample);    
    [y,sorder] = sort(pvals);
    y = -1*log10(y);
else
    step = (max(pvals) - min(pvals))/1000;
    sample = min(pvals):step:max(pvals);
    x = sample;
    [y,sorder] = sort(pvals);
end
clear sample;
figure();
hold on;
grid;
%maxh=ceil(log10(length(x)));
if (docolor)
    colormap('jet');
    cout = colorscale(sorder);
    scatter(x,y,50,cout);
else
    plot(x,y,'bx');
end
plot(x,x,'r-');

% highlight those hits meeting FDR <= 5%
[porder,sorder] = sort(pvals);
fdr = mafdr(porder);
hits = fdr < fdrthresh; % 2
hold on;
scatter(x(hits),y(hits),'ro','SizeData',60);

% write labels:
if (~isempty(pointlabels))
    ps = pointlabels(sorder);
    ps = strcat(ps,'  .');
    ht = text(x(hits),y(hits),ps(hits),'HorizontalAlignment','right');
    set(ht,'Rotation',-35)
end



% compute genomic inflation factor:
gif = (x * y') / (x * x');
xlabel('Expected -log10 p-values');
ylabel('Observed -log10 p-values');
if (~nogif)
    text(max(x)/50, max(y), ['  GIF = ',num2str(gif)], 'FontSize',14);
end
% compute lambda:
lambda = lambdaGIF(pvals);
if (~nogif)
    text(max(x)/50, max(y) - max(y)/12, ['  lambda = ',num2str(lambda)], 'FontSize',14);
end
clear pvals;


