% Script generated by Brainstorm (14-Jun-2021)

% Input files
sFiles = {...
    'link|@default_study/results_MN_EEG_KERNEL_200510_1625.mat|GGroup2-ND0079_rest_postcomp/@rawGGroup2-ND0079_rest_postcomp/data_0raw_GGroup2-ND0079_rest_postcomp.mat'};

% Start a new report
bst_report('Start', sFiles);

% Process: Time-frequency (Morlet wavelets)
sFiles = bst_process('CallProcess', 'process_timefreq', sFiles, [], ...
    'clusters',      {}, ...
    'scoutfunc',     1, ...  % Mean
    'edit',          struct(...
         'Comment',         'Power,FreqBands', ...
         'TimeBands',       [], ...
         'Freqs',           {{'delta', '2, 4', 'mean'; 'theta', '5, 7', 'mean'; 'alpha', '8, 12', 'mean'; 'beta', '15, 29', 'mean'; 'gamma1', '30, 59', 'mean'; 'gamma2', '60, 90', 'mean'}}, ...
         'MorletFc',        1, ...
         'MorletFwhmTc',    3, ...
         'ClusterFuncTime', 'none', ...
         'Measure',         'power', ...
         'Output',          'all', ...
         'SaveKernel',      0), ...
    'normalize2020', 0, ...
    'normalize',     'none');  % None: Save non-standardized time-frequency maps

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);
% bst_report('Export', ReportFile, ExportDir);

