function fx_exportBstTimeFreqPowByTag(type, tag, csvfile, expectedSubjects)

%%== RELATIVE POWER
sFiles = fx_getBstTimeFreqFiles();
expected_total_subjects = expectedSubjects;
sFiles = fx_getBstSelectByTag(sFiles, ...
    tag, 'select', expected_total_subjects);

% export subject level values
switch type
    case 'source'
        process_extract_values_wrapper(sFiles, csvfile, expected_total_subjects);
    case 'electrode'
        [sResultsPow, VariableNames] = ...
            fx_BstExtractValuesElec(sFiles);
        fx_customSaveResultsCSV(csvfile, sResultsPow, VariableNames);
end

end