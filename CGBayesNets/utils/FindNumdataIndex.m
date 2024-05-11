function index = FindNumdataIndex(columnname,cols,numcolindex)
% index = FindNumdataIndex(columnname,cols,numcolindex)
%
% Function to find the index of a column given it's name as a string, from
% the output of an RCSVLoad() function.
%
% INPUT :
%   COLUMNNAME : this is a string, the name of the column for which the
%       index shold be returned.
%   COLS : the cellarray of names of columns, from the output of
%       RCSVLoad()
%   NUMCOLINDEX : either the variable NUMCOLINDEX or STRCOLINDEX returned
%       from the RCVSLoad() function, depending on whether or not COLUMNNAME
%       represents a string or a numeric data type.
%
% OUTPUT : 
%   INDEX : the index into NUMDATA or STRDATA corresponding to COLUMNNAME.
%       NUMDATA and STRDATA are output from RCSVLoad();
%
% (c) Michael McGeachie 2018.  

% find the index:
ind = find(strcmp(columnname,cols));
index = find(numcolindex == ind);

