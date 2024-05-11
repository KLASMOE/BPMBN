%% misc testing:
% common parameter values:
priorPrecision.nu = 25;
priorPrecision.sigma2 = 1;
priorPrecision.alpha = 25;
priorPrecision.maxParents = 3;

%% now learn the chess example using LearnStructure:
[data, cols] = ReadHeaderDataFile('chess-kr-vs-kp.txt');
MBNet = LearnStructure(data,cols,'class',priorPrecision, 'chess-net');
MBNet = LearnParamsBN(MBNet);
disctables = GetProbTables(MBNet);

%%
[data, cols] = ReadHeaderDataFile('winedata.txt');
MBNet = LearnStructure(data, cols, 'class', priorPrecision, 'winetest-net');
MBNet = LearnParamsBN(MBNet);
conttables = GetProbTables(MBNet);

