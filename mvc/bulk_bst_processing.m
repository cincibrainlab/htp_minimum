% 20200828-Source-Manuscript-Figure.m
% 2020-08-f27
% E. Pedapati

% Code Outline
% 1. 
bst_fullfile(DataSrc, '@intra', '*.mat')
%% Calculate ABS power
bst_get('DirAnalysisInter')
bst_get('AnalysisInterStudy')         
% get all TF files
sFiles = bst_process('CallProcess', ...
    'process_select_files_matrix', [], [],...
    'includeintra',  1,...
    'includecommon', 1);
sStudy = bst_get('Study');
bst_get('GetFilenames',1,2,'pdata')
sFiles_stat = bst_process('CallProcess', 'process_select_tag', sFiles, [], ...
'tag',    'sFiles_ABS_Raw_Male', ...
'search', 1, ...
'select', 1); 
% Process: Select files using search query: multiply
sFiles = bst_process('CallProcess', 'process_select_search', [], [], ...
    'search',     'sFiles_', ...
    'includebad', 1);

sFiles = bst_process('CallProcess', ...
    'process_select_files_results', [], []);
sFiles = bst_process('CallProcess', ...
    'process_select_files_data', [], []);
sFiles = bst_process('CallProcess', ...
    'process_select_files_matrix', [], []);
sAllData = bst_process('CallProcess', 'process_select_files_data', [], [], ...
     'includeintra',  1);

DirIntra = bst_get('DirAnalysisIntra');
 sAvgSrc = bst_process('CallProcess', 'process_select_files_data', [], [], ...
        'subjectname',   'All', ...
        'condition', DirIntra, ...
        'includeintra',  1);
    
      sAvgSrc = bst_process('CallProcess', 'process_select_files_data', [], [], ...
        'subjectname',   'All', ...
       'condition',     DirIntra, ...
        'includeintra',  1);

% Process: Select files using search query: multiply
sFiles = bst_process('CallProcess', 'process_select_search', [], [], ...
    'search',     'PSD', ...
    'includebad', 1);


% create specific datasets
sFiles_stat = bst_process('CallProcess', 'process_select_tag', ... [], ...
'tag',    'sFiles_', ...
'search', 1, ...
'select', 1); 

% create specific datasets
sFiles_raw = bst_process('CallProcess', 'process_select_tag', ...
    sFiles, [], ...
'tag',    'PSD: 79/2000ms Power,FreqBands,ABS', ...
'search', 2, ...
'select', 1); 

sFiles_raw2 = bst_process('CallProcess', 'process_select_tag', ...
    sFiles_raw, [], ...
'tag',    'multiply', ...
'search', 2, ...
'select', 2); 

cfg = [];
cfg.fwhm = 2;
cfg.overwrite = 0;

sFiles_smooth = bst_process('CallProcess', ...
    'process_ssmooth_surfstat', sFiles_raw2, [], ...
    'fwhm', cfg.fwhm, ...
    'overwrite', cfg.overwrite);

sFiles_norm = bst_process('CallProcess', 'process_select_tag', ...
    sFiles, [], ...
'tag',    'multiply', ...
'search', 2, ...
'select', 1); 


cfg = [];
cfg.timewindow = [0, 80];
cfg.win_length = 2;
cfg.win_overlap = 50;
cfg.scoutfunc = 1;
cfg.freqs = {'delta', '2, 4', 'mean'; 'theta', '4, 7.5', 'mean'; 'alpha1', '8, 10', 'mean'; 'alpha2', '10.5, 12.5', 'mean'; 'beta', '15, 29', 'mean'; 'gamma_lo', '30, 55', 'mean'; 'gamma_hi', '63, 90', 'mean'};

sFiles = bst_process('CallProcess', 'process_psd', sFiles, [], ...
    'timewindow', cfg.timewindow, ...
    'win_length', cfg.win_length, ...
    'win_overlap', cfg.win_overlap, ...
    'clusters', {}, ...
    'scoutfunc', cfg.scoutfunc, ...% Mean
    'win_std', 0, ...
    'edit', struct(...
    'Comment', 'Power,FreqBands,ABS', ...
    'TimeBands', [], ...
    'Freqs', {cfg.freqs}, ...
    'ClusterFuncTime', 'none', ...
    'Measure', 'power', ...
    'Output', 'all', ...
    'SaveKernel', 0));

cfg = [];
cfg.fwhm = 3;
cfg.overwrite = 0;
sFile = o.bst_spatialSmoothing(cfg, sFile);
tag = o.bst_generateTag(r, sprintf('ssmooth_%s', num2str(cfg.fwhm)));

sFiles = bst_process('CallProcess', 'process_ssmooth_surfstat', sFiles, [], ...
    'fwhm', cfg.fwhm, ...
    'overwrite', cfg.overwrite);

% Process: Add tag: test
sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
    'tag', 'abs7326', ...
    'output', 1);

%%
% Process: Spectrum normalization
sFiles_norm = bst_process('CallProcess', 'process_tf_norm', sFiles, [], ...
    'normalize', 'multiply2020', ...  % 1/f compensation (multiply power by frequency)
    'overwrite', 1);


 %% create Grp Ids
grpIdxTmp = {sFiles.Condition}';
for grpi = 1 : length(grpIdxTmp)
    if contains(grpIdxTmp{grpi}, 'rawGGroup2')
        grpIdx(grpi) = 0;
    else
        grpIdx(grpi) = 1;
    end
end

sFiles(grpIdx ==1)
sFiles(grpIdx ==0)

%% create Grp Ids with Sex
clindata = readtable("C:\Users\ernie\Dropbox\NN_RestingSource\Results\data-raw\fxssource_table_01_demographics_data.csv")

eegids = clindata.eegid;
sids = clindata.sex;
gids = clindata.group;
clinidx = [];
for clini = 1 : length(eegids)
    
    tmpidx = find(contains( grpIdxTmp, eegids{clini}));
    
    if isnumeric(tmpidx) & ~isempty(tmpidx)
        clinidx{tmpidx} = [sids{clini} '_' gids{clini}];
    end
    
end

for grpi = 1 : length(clinidx)
    
    switch clinidx{grpi}
        case 'M_FXS'
          grpIdxClin(grpi) = 1;
        case 'M_TDC'
            grpIdxClin(grpi) = 2;
        case 'F_FXS'
            grpIdxClin(grpi) = 3;
        case 'F_TDC'
            grpIdxClin(grpi) = 4;
    end
end



%% Source Stats - ALL SUBJECTS
cfg = [];
%cfg.r = am.bst_generateRandomTag;
cfg.tag = 'abs_src_stats';
cfg.timewindow       = [0 80];
cfg.isabs            = 0;
cfg.avgtime          = 1;
cfg.avgfreq          = 0;
cfg.randomizations   = 2000;
cfg.statisticstype   = 1;
cfg.tail             = 'two';
cfg.correctiontype   = 1;  % 1 default 2 cluster
cfg.minnbchan        = 1;
cfg.clusteralpha     = 0.05;

sFiles = sFiles_smooth

cfg.tagstr = 'sFiles_ABS_Raw_Male';
sFiles_ABS_Raw_Male = ...
    process_ft_sourcestatistics_wrapper(...
    sFiles(grpIdxClin == 1), ...
    sFiles(grpIdxClin == 2), ...
    cfg);

cfg.tagstr = 'sFiles_ABS_Raw_Female';
sFiles_ABS_Raw_Female = ...
    process_ft_sourcestatistics_wrapper(...
    sFiles(grpIdxClin == 3), ...
    sFiles(grpIdxClin == 4), ...
    cfg);

cfg.tagstr = 'sFiles_ABS_Raw_All';
sFiles_ABS_Raw_All = ...
    process_ft_sourcestatistics_wrapper(...
    sFiles(grpIdxClin == 1 | grpIdxClin == 3), ...
    sFiles(grpIdxClin == 2 | grpIdxClin == 4), ...
    cfg);

cfg.tagstr = 'sFiles_ABS_Norm_Male';
sFiles_ABS_Norm_Male = ...
    process_ft_sourcestatistics_wrapper(...
    sFiles_norm(grpIdxClin == 1), ...
    sFiles_norm(grpIdxClin == 2), ...
    cfg);

cfg.tagstr = 'sFiles_ABS_Norm_Female';
sFiles_ABS_Norm_Female = ...
    process_ft_sourcestatistics_wrapper(...
    sFiles_norm(grpIdxClin == 3), ...
    sFiles_norm(grpIdxClin == 4), ...
    cfg);

cfg.tagstr = 'sFiles_ABS_Norm_All';
sFiles_ABS_Norm_All = ...
    process_ft_sourcestatistics_wrapper(...
    sFiles_norm(grpIdxClin == 1 | grpIdxClin == 3), ...
    sFiles_norm(grpIdxClin == 2 | grpIdxClin == 4), ...
    cfg);

% get p-values for each group
statResultFiles = {sFiles_ABS_Norm_Male, ...
sFiles_ABS_Norm_Female, ...
sFiles_ABS_Norm_All, ...
sFiles_ABS_Raw_Male, ...
sFiles_ABS_Raw_Female, ...
sFiles_ABS_Raw_All};

statResultsCsv = {};
for stati = 1 : length(statResultFiles)
    
    statResultsCsv{stati, 1} = statResultFiles{stati}.Comment;
    statResultsCsv{stati, 2} = 'fdr';
    statResultsCsv{stati, 3} = 0.05;
    statcsd =  in_bst_data(statResultFiles{stati}.FileName)
    
    % get p-map and define treshold p for multiple comparisons
    pmap = statcsd.pmap;
    StatThreshOptions.pThreshold =statResultsCsv{stati, 3};
    StatThreshOptions.Correction =  statResultsCsv{stati, 2};
    StatThreshOptions.Control = [1 2 3];
    [pmask, pthresh] = bst_stat_thresh(pmap(:,1,:), StatThreshOptions);
        statResultsCsv{stati, 4} = pthresh;

end


%%           
process_extract_values_wrapper(sFiles, '210610_ABS_P1_ALL.csv', 141);
process_extract_values_wrapper(sFiles_norm, '210610_ABS_P1_ALL_norm.csv', 141);

%%
count = 1;
TF = {};
for stati = 1 : length(statResultFiles)
    currentstat = statResultFiles{stati};
    
   % TF = cell(length(sFiles_norm), 9);
    for sigi = 1 : length(sFiles)
        eegid = extractBetween(sFiles_norm(sigi).SubjectName, '-N','_p');
        groupid = extractBetween(sFiles_norm(sigi).SubjectName, 'G','-N');
        TF{count, 1} = eegid{1};
        TF{count, 2} = groupid{1};
        
        rawvalue = in_bst_data(sFiles_norm(sigi).FileName);
        
        TF{count, 3} = currentstat.Comment;
        TF{count, 4} = rawvalue.Method;
        TF{count, 5} = 'statsig';
        
            statcsd =  in_bst_data(statResultFiles{stati}.FileName);
    
            % get p-map and define treshold p for multiple comparisons
            pmap = statcsd.pmap;
            tmap = squeeze(statcsd.tmap(:,1,:));
            pos_tmap = tmap > 0;
            neg_tmap = tmap < 0;
            StatThreshOptions.pThreshold =  statResultsCsv{stati, 3};
            StatThreshOptions.Correction =  statResultsCsv{stati, 2};
            StatThreshOptions.Control = [1 2 3];
            [pmask, pthresh] = bst_stat_thresh(pmap(:,1,:), StatThreshOptions);
            %statResultsCsv{stati, 4} = pthresh;
            
            sigvalues = squeeze(rawvalue.TF) .* squeeze(pmask);
            possigvalues = sigvalues .* pos_tmap;
            negsigvalues = sigvalues .* neg_tmap;
            possigvalues(possigvalues == 0) = NaN;
            negsigvalues(negsigvalues == 0) = NaN;

            TF(count, 6:12) = num2cell(mean(possigvalues, 1, 'omitnan'));
            TF(count, 13:19) = num2cell(mean(negsigvalues, 1, 'omitnan'));
            count = count + 1;
    end
    
end
%%
% statresfiles
% band
% label
% pos
% neg
% counts

% needs frequency Freqs cell array

count = 1;
roi_res = {};
for stati = 1 : length(statResultFiles)
    
    currentstat = statResultFiles{stati};
    statcsd =  in_bst_data(currentstat.FileName);
    
    roistrip = 1 : length(atlas.Scouts);
    roilabels = {atlas.Scouts.Label};
    roiregions = {atlas.Scouts.Region};
    
    pmap = squeeze(statcsd.pmap(:,1,:));
    tmap = squeeze(statcsd.tmap(:,1,:));
    pos_tmap = tmap > 0;
    neg_tmap = tmap < 0;
    
    StatThreshOptions.pThreshold = 0.05;
    StatThreshOptions.Correction = 'fdr';
    StatThreshOptions.Control = [1 2 3];
    
    [pmask, pthresh] = bst_stat_thresh(pmap, StatThreshOptions);
    
    sigvalues = squeeze(rawvalue.TF) .* squeeze(pmask); % only include sig values
    possigvalues = sigvalues .* pos_tmap; % get only positive direction
    negsigvalues = sigvalues .* neg_tmap; % get only negative direction
    possigvalues(possigvalues == 0) = NaN; % important for mean function
    negsigvalues(negsigvalues == 0) = NaN;
    
    for freqi =  1 : length(Freqs)
        
        siglabelspos = find(~isnan(possigvalues(:,freqi)));
        siglabelsneg = find(~isnan(negsigvalues(:,freqi)));
        
        for vert_i = 1 : length(siglabelspos)
            if ~isempty(siglabelspos
            roi_res_pos{vert_i} = vertex2roi(atlas, siglabelspos(vert_i));
            roi_res_neg{vert_i} = vertex2roi(atlas, siglabelsneg(vert_i));
        end
        

        
        labelcnts= categorical(cell2mat(roi_res_neg),roistrip,roilabels)
        regioncnts = categorical(cell2mat(roi_res_neg),roistrip,roiregions)
        
        [roilabels', num2cell(countcats(labelcnts)')]
        [roilabels', num2cell(countcats(regioncnts)')]

        
        for labeli = 1 : length(atlas.Scouts)
            sel_scout = atlas.Scouts(labeli);
            
            

                 
            roi_res{count, 1} = currentstat.Comment;
            roi_res{count, 2} = Freqs{freqi, 1};
            roi_res{count, 3} = sel_scout.Label;
            count = count + 1;
            
        end
    end
    
end
roi_res(1:3,:)
%%

for roi_i = 1 : length(Freqs)
    
    roi_search_idx = find(neg_tmap(:,4)==1);
    for vert_i = 1 : length(roi_search_idx)
        roi_res{vert_i} = vertex2roi(atlas, roi_search_idx(vert_i));
    end
    
end

roistrip = 1 : numel(scouts);
roilabels = {scouts.Label};
roiregions = {scouts.Region};

labelcnts= categorical(cell2mat(roi_res),roistrip,roilabels)
regioncnts = categorical(cell2mat(roi_res),roistrip,roiregions)

[roilabels', num2cell(countcats(labelcnts)')]
[roilabels', num2cell(countcats(regioncnts)')]

summary(labelcnts)
summary(regioncnts)

roicats = categorical(cell2mat(roi_res));
categories(roicats)
U = unique(cateogrical(cell2mat(roi_res)));
N = histc(roi_res,U);
[U,N]
%%


  fun = @(x) scouts(x).Vertices == 15002 % useful for complicated fields
  tf2 = arrayfun(fun, 1:numel(atlas),'uni',0)
  index2 = find(tf2)

Freq = [{'eegid','group','Compare','Method','Type'} rawvalue.Freqs(:,1)'];
sig_result_table = cell2table(TF, 'VariableNames', Freq);

%%


% get time frequency values variable
subcsd = in_bst(sFiles_values.FileName); % Freqs: {8×3 cell} TF: [15002×141×8 double] Time: [1×141 double]
sCortex = in_tess_bst(subcsd.SurfaceFile);
subcsd.iAtlas = find(strcmpi({sCortex.Atlas.Name}, 'Desikan-Killiany'));
subcsd.Atlas = sCortex.Atlas(subcsd.iAtlas)

%%


%% Create relative power spreadsheet from each source map

% load protocol P1_FXS70_71_NOBZO (via GUI)
% brainstorm;

% gather source files (TF analysis should be complete)
% make sure assigned tag to analysis is unique, i.e. no_bzo_relpow_epsilon 
sFiles = bst_process('CallProcess', 'process_select_files_results', [], []);
 
% Gather all time frequency files for each subject based on tag
% Process: Select timefreq files in: */*/no_bzo_relpow_epsilon 
sFiles = bst_process('CallProcess', 'process_select_files_timefreq', sFiles, [], ...
    'subjectname',   '', ...
    'condition',     '', ...
    'tag',           'ssmooth3', ...
    'includebad',    0, ...
    'includeintra',  0, ...
    'includecommon', 0);

% extract all values at each vertex
% Process: Extract values: [-1.000s,80.000s] 2-140Hz
sFiles = bst_process('CallProcess', 'process_extract_values', sFiles, [], ...
    'timewindow', [0, 80], ...
    'freqrange',  [2, 140], ...
    'rows',       '', ...
    'isabs',      0, ...
    'avgtime',    1, ...
    'avgrow',     0, ...
    'avgfreq',    0, ...
    'matchrows',  0, ...
    'dim',        2, ...  % Concatenate time (dimension 2)
    'Comment',    '');

% get time frequency values variable
subcsd = in_bst(sFiles.FileName); % Freqs: {8×3 cell} TF: [15002×141×8 double] Time: [1×141 double]
sCortex = in_tess_bst(subcsd.SurfaceFile);
subcsd.iAtlas = find(strcmpi({sCortex.Atlas.Name}, 'Desikan-Killiany'));
subcsd.Atlas = sCortex.Atlas(subcsd.iAtlas);

specialScout = subcsd.Atlas.Scouts; 
atlas = subcsd.Atlas;  
specialScout= specialScout(1);  % make it single region
template = specialScout(end);
sigfdr =  0.0211;
scoutindex = length(specialScout);

for i = 1 : size(statcsd.pmap,3) % 8 bands of sig
    verts = squeeze(statcsd.pmap(:,1,i)) < sigfdr; % index of sig    
    freqlabel = subcsd.Freqs{i,1};  % get label of band
    scoutindex = i;
    template.Vertices = find(verts);  % get index
    template.Label = freqlabel;
    template.Region = 'StatSig';
    specialScout(scoutindex) = template;
end
%%
%output table format
% ID EEGID GROUP BAND SIGPOWER
raw_array_of_subj_names = subcsd.History(5:end,3);
expected_total_subjects = 141;  % sanity check for correct array

% preallocation
primaryid = 1;  % primary table key
result_cell = {expected_total_subjects*size(specialScout,2)};

% check # of subjects
assert(expected_total_subjects == length(raw_array_of_subj_names))

for i = 1 : expected_total_subjects
       
       % extract group and subject name
       slabel_tmp = strsplit(raw_array_of_subj_names{i});
       slabel_tmp2 = slabel_tmp(end);
       slabel_tmp3 = strsplit(slabel_tmp2{1},'/');
       slabel_tmp4 = slabel_tmp3{1};
       group_eegid = strsplit(slabel_tmp4,'-');
       eegid = erase(group_eegid{2}, '_postcomp');
       eegid = regexprep(eegid, 'N','');
       groupid = group_eegid{1};
       
       for j = 1 : size(specialScout,2)
           
           sig_vertices = specialScout(j).Vertices;
           bandlabel = specialScout(j).Label;
           region = specialScout(j).Region;
           avgTF = mean(subcsd.TF(sig_vertices,i,j));

           result_cell{primaryid,1} = primaryid;
           result_cell{primaryid,2} = groupid;
           result_cell{primaryid,3} = eegid;
           result_cell{primaryid,4} = bandlabel;
           result_cell{primaryid,5} = region;
           result_cell{primaryid,6} = avgTF;
           result_cell{primaryid,7} = bandlabel;
           
           primaryid = primaryid + 1;
       end
end                                                                                                                                                                                                                             

size(result_cell)

 result_table = cell2table(result_cell,'VariableNames', {'id', 'group', ...
       'eegid','label','region','value','bandname'});

   writetable(result_table, ...
       fullfile(p.am.sub(1).pathdb.analysis,'20200906_relpow_sigonly_source_noepsilon.xlsx'));
   
%%
% all subjects
results = {};
subjectnames = subcsd.History(5:end,3);


count = 1;
if length(subjectnames) == 141 % hard coded subject number for sanity check
    for i = 1 : length(subjectnames)
        subname = subjectnames{i};
        % fix name
            subnametmp = strsplit(subname);
            subnametmp2 = subnametmp(end);
            tmpstr = strsplit(subnametmp2{1},'/');
            tmpstr = tmpstr{1};            
            bst_file = strsplit(tmpstr,'-');
            eegid = erase(bst_file{2}, '_postcomp');
            eegid = regexprep(eegid, 'N','');
            groupid = bst_file{1};


        for j = 1 : size( atlas.Scouts,2)
            
            sname =  atlas.Scouts(j).Label;

            sregion =  atlas.Scouts(j).Region;
            svertex =  atlas.Scouts(j).Vertices;
            Freqs = subcsd.Freqs;
            for k = 1 : size(Freqs,1)
                avgTF = mean(subcsd.TF(svertex,i,k));
                bandname = Freqs{k,1};
                results{count,1} = count;
                results{count,2} = groupid;
                results{count,3} = eegid;
                results{count,4} = sname;
                results{count,5} = sregion;
                results{count,6} = bandname;
                results{count,7} = avgTF;
                
                count = count + 1;
            end
        end
        
    end
end


bstsourcepow = cell2table(results, 'VariableNames', {'id','group',...
    'eegid','label','region','bandname','value'});

writetable(bstsourcepow, fullfile(p.am.sub(1).pathdb.analysis,'20200906_relpow_source_noepsilon.xlsx'));

%% create spectrogram from source data
% legacy not used
% get variable "noband" from GUI FFT 2 to 250 Hz

plot(squeeze(mean(mean(noband.TF,1),2)))
%%
count =1
results = {};
subjectnames = noband.History(5:end,3);

if length(subjectnames) == 141 % hard coded subject number for sanity check
    for i = 1 : length(subjectnames)
        subname = subjectnames{i};
        % fix name
        subnametmp = strsplit(subname);
        subnametmp2 = subnametmp(end);
        tmpstr = strsplit(subnametmp2{1},'/');
        tmpstr = tmpstr{1};            
        bst_file = strsplit(tmpstr,'-');
        eegid = erase(bst_file{2}, '_postcomp');
        eegid = regexprep(eegid, 'N','');
        groupid = bst_file{1};
        
          for j = 1 : size( atlas.Scouts,2)
            
            sname =  atlas.Scouts(j).Label;

            sregion =  atlas.Scouts(j).Region;
            svertex =  atlas.Scouts(j).Vertices;
            Freqs = noband.Freqs;
            for k = 1 : length(Freqs)
                avgTF = mean(noband.TF(svertex,i,k));
                bandname = Freqs(k);
                results{count,1} = count;
                results{count,2} = groupid;
                results{count,3} = eegid;
                results{count,4} = sname;
                results{count,5} = sregion;
                results{count,6} = bandname;
                results{count,7} = avgTF;
                
                count = count + 1;
            end
        end
        
        
    end
    
end

bstsourcepow_noband = cell2table(results, 'VariableNames', {'id','group',...
    'eegid','label','region','bandname','value'});

writetable(bstsourcepow_noband, fullfile(p.am.sub(1).pathdb.analysis,'BSTPOWER082020_noband.csv'));

    