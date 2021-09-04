function resultTable = fx_customStatReportSigLabels( sStats, atlas )

roi_res = {};

roistrip = 1 : length(atlas.Scouts);
roilabels = {atlas.Scouts.Label};
roiregions = {atlas.Scouts.Region};

StatThreshOptions.pThreshold = 0.05;
StatThreshOptions.Correction = 'fdr';
StatThreshOptions.Control = [1 2 3];

runningTable = cell(1,5);
runningTable_tmp ={};

for stati = 1 : length(sStats)
    
    currentStat = sStats(stati);
    Freqs = currentStat.Freqs;
    currentStatComment = currentStat.Comment;
    
    pmap = squeeze(currentStat.pmap(:,1,:));
    tmap = squeeze(currentStat.tmap(:,1,:));
    pos_tmap = tmap > 0;
    neg_tmap = tmap < 0;
    
    [pmask, pthresh] = bst_stat_thresh(pmap, StatThreshOptions);
    
    sigvalues = squeeze(pmask); % only include sig values
    possigvalues = sigvalues .* pos_tmap; % get only positive direction
    negsigvalues = sigvalues .* neg_tmap; % get only negative direction
    
    possigvalues(possigvalues == 0) = NaN; % important for mean function
    negsigvalues(negsigvalues == 0) = NaN;
    
    for freqi =  1 : length(Freqs)
        currentFreq = Freqs(freqi,1);
        roi_res_pos = {};
        roi_res_neg = {};
        
        siglabelspos = find(~isnan(possigvalues(:,freqi)));
        siglabelsneg = find(~isnan(negsigvalues(:,freqi)));
        
        if ~isempty(siglabelspos)
            for vert_i = 1 : length(siglabelspos)
                roi_res_pos{vert_i} = vertex2roi(atlas, siglabelspos(vert_i));
            end
        else
            roi_res_pos = num2cell(repmat(0, 1, length(pmap)));
        end
        labelcats_pos = categorical(cell2mat(roi_res_pos),roistrip,roilabels);
        labelcnts_pos = [roilabels', num2cell(countcats(labelcats_pos)')];
        
        if ~isempty(siglabelsneg)
            for vert_i = 1 : length(siglabelsneg)
                roi_res_neg{vert_i} = vertex2roi(atlas, siglabelsneg(vert_i));
            end
        else
            roi_res_neg = num2cell(repmat(0, 1, length(pmap)));
        end
        labelcats_neg = categorical(cell2mat(roi_res_neg),roistrip,roilabels);
        labelcnts_neg = [roilabels', num2cell(countcats(labelcats_neg)')];
        
        runningTable_tmp = [repmat({currentStatComment},length(roilabels'),1) repmat(currentFreq,length(roilabels'),1) labelcnts_pos labelcnts_neg(:,2)];
        runningTable = [runningTable; runningTable_tmp];
    end
    
    
end

resultTable = runningTable(2:end,:);
end