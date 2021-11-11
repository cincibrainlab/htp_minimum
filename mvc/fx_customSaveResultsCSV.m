function res = fx_customSaveResultsCSV( outfile_csv, results_cell, VariableNames )

bstsourcepow = cell2table(results_cell, 'VariableNames', VariableNames);

writetable(bstsourcepow, fullfile(outfile_csv));
res = bstsourcepow;
end