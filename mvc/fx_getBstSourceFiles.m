function sFiles = fx_getBstTimeFreqFiles()
sFiles = bst_process('CallProcess', ...
    'process_select_files_results', [], [],...
    'includeintra',  1,...
    'includecommon', 1);
str = sprintf('Number of files selected: %d', numel(sFiles));
disp(str);
end