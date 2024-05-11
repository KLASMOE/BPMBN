% Script: CGBayesNets are used to load the known network structure, and output the computation time.

%% CGBayesNets setting
run('./bnpathscript');

% common parameter values:
priorPrecision.nu = 1;
priorPrecision.sigma2 = 1;
priorPrecision.alpha = 10;
priorPrecision.maxParents = 3;

%% Simulation setting
analysis_title = 'CompareBPMBN';
pheno = 'V1';

nodes = [50];
continuous = [0, 0.25, 0.5, 0.75, 1];
samples = [1];
average_times = zeros(length(nodes), length(continuous), 10, length(samples));

n_networks = 100;
totalIterations = length(nodes) * length(continuous) * 10 * length(samples) * n_networks;
currentIteration = 0;

for n_idx = 1:length(nodes)
    n_node = nodes(n_idx);
    for p_idx = 1:length(continuous)
        p_continuous = continuous(p_idx);
        for e = 0:9
            n_evidence = floor((e * n_node) / 10);
            for s_idx = 1:length(samples)
                n_sample = samples(s_idx);

                times = zeros(1, n_networks);
                fprintf(1, 'Computing nodes=%s continuous=%s evidence=%s samples=%s\n', ...
                    num2str(n_node), num2str(p_continuous), num2str(n_evidence), num2str(n_sample));
                for i = 1:n_networks
                    % Construct filename
                    base_dir = '../BPMBN/BPMBN/results';
                    formatted_continuous = num2str(p_continuous*100, '%.0f'); % Formats the continuous value
                    data_filename = sprintf('%s/data/data_V%dC%sE%dN%d_%d.csv', base_dir, n_node, formatted_continuous, n_evidence, n_sample, i-1);
                    stru_filename = sprintf('%s/structure/network_V%dC%s.m', base_dir, n_node, formatted_continuous);
                    prob_filename = sprintf('%s/structure/structure_V%dC%s_%d.txt', base_dir, n_node, formatted_continuous, i-1);
                    evid_filename = sprintf('%s/evidence/evidence_V%dC%sE%d.csv', base_dir, n_node, formatted_continuous, n_evidence);

                    % Load data
                    [data, cols] = RCSVLoad(data_filename, false);

                    % Load network
                    run(stru_filename);
                    nodesketch = our;

                    % Learn params
                    ourBN = BayesNet(nodesketch, analysis_title, [],[],[],[],[],data,cols, pheno,priorPrecision);
                    ourBN = LearnParamsBN(ourBN);

                    % Check if there are any evidence variables
                    if n_evidence == 0
                        evinds = false(1, length(ourBN.cols));
                    else
                        [evi_data, evi_cols] = RCSVLoad(evid_filename, false);
                        evidenceIndices = cellfun(@str2double, evi_cols);
                        evinds = false(1, length(ourBN.cols));
                        evinds(evidenceIndices) = true;
                    end
                    
                    % Computation
                    predictionIndices = find(~evinds);
                    tic; 
                    for j = predictionIndices
                        predvars = j;
                        evvars = find(evinds);
                        [preds, z, lep] = PredictPartial(ourBN, evvars, predvars);
                    end
                    times(i) = toc;
                    disp(times(i));
                    currentIteration = currentIteration + 1;
                    fprintf('Progress：%.2f%%\n', (currentIteration / totalIterations) * 100);
                end
                
                % Save all results
                average_times(n_idx, p_idx, e+1, s_idx) = mean(times);
                save(fullfile('V50N1B100.mat'), 'average_times'); 
            end
        end
    end
end


