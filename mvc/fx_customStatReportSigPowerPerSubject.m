function resultTable = fx_customStatReportSigPowerPerSubject( p, sStats, sFiles )
[subnames, groupids] = fx_customGetSubNames(  sFiles, p, 'default' );

roi_res = {};

%roistrip = 1 : length(atlas.Scouts);
%roilabels = {atlas.Scouts.Label};
%roiregions = {atlas.Scouts.Region};

StatThreshOptions.pThreshold = 0.05;
StatThreshOptions.Correction = 'fdr';
StatThreshOptions.Control = [1 2 3];

runningTable = cell(1,6);
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
        
        
        % load subject data arrays
        for subi = 1 : length(sFiles)
            if ~isempty(siglabelspos)
                currentSub = sFiles(subi);
                subdata = in_bst(currentSub.FileName);
                currentTF = squeeze(subdata.TF(:,1,freqi));
                meanPowAtSigPosIndex{subi} = mean(currentTF(siglabelspos));
                
            else
                meanPowAtSigPosIndex{subi} = NaN;
            end
            
        end
        
        
        % load subject data arrays
        for subi = 1 : length(sFiles)
            if ~isempty(siglabelsneg)
                currentSub = sFiles(subi);
                eegid{subi} = subnames(subi);
                groupid{subi} = groupids(subi);
                subdata = in_bst(currentSub.FileName);
                currentTF = squeeze(subdata.TF(:,1,freqi));
                meanPowAtSigNegIndex{subi} = mean(currentTF(siglabelsneg));
            else
                meanPowAtSigNegIndex{subi} = NaN;
                
            end
        end
        
        tableLength = length(sFiles);
        runningTable_tmp = [ ...
            repmat({currentStatComment},tableLength,1) ...
            subnames' ...
            groupids' ...
            repmat(currentFreq,tableLength,1) ...
            meanPowAtSigPosIndex' ...
            meanPowAtSigNegIndex'];
        runningTable = [runningTable; runningTable_tmp];
    end
end

resultTable = runningTable(2:end,:);
end