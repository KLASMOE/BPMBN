function newcpt = MultiFactorSum(cpt, varname)
% 
%
% use factorsum to get rid of all the variables in the CPT EXCEPT the
% varname.


% loop through all vars in the CPT
newcpt = cpt;
done = false;
while (~done)
    done = true;
    % find a factor to reduce:
    for i = 1:length(newcpt.factors)
        % if they're not VARNAME and not already reduced to a single value, call
        % factorsum
        if (~strcmpi(newcpt.factors{i}, varname) && length(newcpt.values{i}) > 1)

            newcpt = factorsum(newcpt,newcpt.factors{i});
            done = false; % check for more factors after this one
            break; % just do one factor at a time
        end
    end
    %check if there are any other factors left:
end