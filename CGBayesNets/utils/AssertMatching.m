function match = AssertMatching(list1, list2)
% function to check two lists for matching elements in order

match = true;

if (iscell(list1))
    DOCELL = true;
else
    DOCELL = false;
end

if (length(list1) ~= length(list2))
    match = false;
    return;
end

for i = 1:length(list1)
    if (DOCELL)
        if (~strcmp(list1{i},list2{i}))
            match = false;
            return;
        end
    else
        if (list1(i) ~= list2(i))
            match = false;
            return;
        end
    end
end

    
    