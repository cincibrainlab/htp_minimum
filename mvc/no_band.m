% Script generated by Brainstorm (14-Jun-2021)

% Input files
sFiles = {...
    'link|@default_study/results_MN_EEG_KERNEL_200510_1625.mat|GGroup2-ND0079_rest_postcomp/@rawGGroup2-ND0079_rest_postcomp/data_0raw_GGroup2-ND0079_rest_postcomp.mat'};

% Start a new report
bst_report('Start', sFiles);

% Process: Power spectrum density (Welch)
sFiles = bst_process('CallProcess', 'process_psd', sFiles, [], ...
    'timewindow',  [0, 80], ...
    'win_length',  2, ...
    'win_overlap', 50, ...
    'units',       'physical', ...  % Physical: U2/Hz
    'clusters',    {'Desikan-Killiany', {'bankssts L', 'bankssts R', 'caudalanteriorcingulate L', 'caudalanteriorcingulate R', 'caudalmiddlefrontal L', 'caudalmiddlefrontal R', 'cuneus L', 'cuneus R', 'entorhinal L', 'entorhinal R', 'frontalpole L', 'frontalpole R', 'fusiform L', 'fusiform R', 'inferiorparietal L', 'inferiorparietal R', 'inferiortemporal L', 'inferiortemporal R', 'insula L', 'insula R', 'isthmuscingulate L', 'isthmuscingulate R', 'lateraloccipital L', 'lateraloccipital R', 'lateralorbitofrontal L', 'lateralorbitofrontal R', 'lingual L', 'lingual R', 'medialorbitofrontal L', 'medialorbitofrontal R', 'middletemporal L', 'middletemporal R', 'paracentral L', 'paracentral R', 'parahippocampal L', 'parahippocampal R', 'parsopercularis L', 'parsopercularis R', 'parsorbitalis L', 'parsorbitalis R', 'parstriangularis L', 'parstriangularis R', 'pericalcarine L', 'pericalcarine R', 'postcentral L', 'postcentral R', 'posteriorcingulate L', 'posteriorcingulate R', 'precentral L', 'precentral R', 'precuneus L', 'precuneus R', 'rostralanteriorcingulate L', 'rostralanteriorcingulate R', 'rostralmiddlefrontal L', 'rostralmiddlefrontal R', 'superiorfrontal L', 'superiorfrontal R', 'superiorparietal L', 'superiorparietal R', 'superiortemporal L', 'superiortemporal R', 'supramarginal L', 'supramarginal R', 'temporalpole L', 'temporalpole R', 'transversetemporal L', 'transversetemporal R'}}, ...
    'scoutfunc',   1, ...  % Mean
    'win_std',     0, ...
    'edit',        struct(...
         'Comment',         'Scouts,Power', ...
         'TimeBands',       [], ...
         'Freqs',           [], ...
         'ClusterFuncTime', 'after', ...
         'Measure',         'power', ...
         'Output',          'all', ...
         'SaveKernel',      0));

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);
% bst_report('Export', ReportFile, ExportDir);
