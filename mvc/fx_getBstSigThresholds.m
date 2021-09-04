function statResultsCsv = fx_getBstSigThresholds(sStats)
statResultsCsv = {};

StatThreshOptions.pThreshold = 0.05;
StatThreshOptions.Correction = 'fdr';
StatThreshOptions.Control = [1 2 3];

for stati = 1 : length(sStats)
    
    current_sStat = sStats(stati);
    
    statResultsCsv{stati, 1} = current_sStat.Comment;
    statResultsCsv{stati, 2} = StatThreshOptions.Correction;
    statResultsCsv{stati, 3} = StatThreshOptions.pThreshold;
    
    [pmask, pthresh] = bst_stat_thresh(current_sStat.pmap(:,1,:), ...
        StatThreshOptions);
    statResultsCsv{stati, 4} = pthresh;
end

end