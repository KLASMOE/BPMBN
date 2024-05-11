function index = getRCSVcolindex(str, cols, numcolindex)
%
% utility function to help with using RCSVLoad() indexing scheme

match = find(strcmp(str, cols));
index = find(match == numcolindex);





