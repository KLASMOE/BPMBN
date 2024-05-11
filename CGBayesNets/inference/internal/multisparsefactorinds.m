function fullinds = multisparsefactorinds(f, varind, varval)
%
% rather than holding a whole table of combinatorial index comutations :
%   (ie, it looks like : 
%       0 0 0 : x1
%       0 0 1 : x2
%       0 1 0 : x3, ... etc )
% this computes the appripriate indices for variables that are assigned
% values varval.
%
% INPUT : 
%   F : CondProbTable
%   VARIND : indices into F.factors{VARIND}
%   VARVAL : value of f.factors{varind} to match VARVAL
%
% returns an index vector into f.probs for values of f.factors{varind}
% matching VARVAL
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.


vlengths = zeros(size(f.values));
for i = 1:length(f.values)
    vlengths(i) = length(f.values{i});
end
cp = cumprod(vlengths);

vrank = zeros(1,length(varind));
for j = 1:length(varind)
    for i = 1:length(f.values{varind(j)})
        if (f.values{varind(j)}(i) == varval(j))
            vrank(j) = i;
            break;
        end
    end
end


fullinds = true(1, cp(end));
for j = 1:length(vrank)
    if (varind(j) > 1)
        fbstart = (vrank(j)-1) * cp(varind(j)-1) + 1;
        fbend = vrank(j) * cp(varind(j)-1);
        blockend = cp(varind(j));
        inds = false(1,blockend);
        inds(fbstart:fbend) = true;

        % slickness: 
        repeats = cp(end) / cp(varind(j));
        inds = repmat(inds,1,repeats);
    else
        inds = false(1,vlengths(varind(j)));
        inds(vrank(j)) = true;
        repeats = cp(end) / cp(varind(j));
        inds = repmat(inds,1,repeats);
    end
    fullinds = fullinds & inds;
end    
    
    
    
    
    
    
    
    
    
    
    
    