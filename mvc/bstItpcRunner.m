%% Optimized ITPC

%% Chirp High Resolution Cortical Investigation
brainstorm;
sStudy = bst_get('Study');
sProtocol = bst_get('ProtocolInfo');
sSubjects =  bst_get('ProtocolSubjects');
sStudyList = bst_get('ProtocolStudies'); 
atlas = fx_BstGetDKAtlasFromSurfaceMat;
sCortex = in_tess_bst('@default_subject/tess_cortex_pial_low.mat');

% get pointers to raw recordings
sFilesSource2 = fx_getBstSourceFiles();
% sFilesSource2 = fx_getBstSelectByTag(sFilesSource, ...
%     'sensorychirp', 'deselect', 75);
% sFilesSource3 = fx_getBstSelectByTag(sFilesSource, ...
%     'MN: EEG(Constr) 2018', 'select', 15505);
% sFilesSource4 = fx_getBstSelectByTag(sFilesSource3, ...
%     '_avgtime', 'deselect', 10850);

[subnames groupids eegid_base] = fx_customGetSubNames(sFilesSource2);
%%
itpc = zeros(15002, length(sFilesSource2));
for i =1 : length(sFilesSource2)
    bstsub = in_bst(sFilesSource2(i).FileName);
    itpc(:,i) = fx_getItpc( bstsub );
end
    % export subject level values
            subcsd.Freqs = {'chirp', '37,43', 'itpc'};
            
            sFiles = fx_getBstTimeFreqFiles();
expected_total_subjects = 75;
sFiles = fx_getBstSelectByTag(sFiles, ...
        'TEMPLATE', 'select', expected_total_subjects);

    csvfile = fullfile('FXSREST_mneItpc.csv');
    process_extract_values_wrapper(  sPow.mneAbspow , csvfile, 75);
  
    groupid_table = innerjoin( cell2table(eegid_base', 'VariableNames', {'eegid'}), clindata(:, {'eegid', 'sex', 'group'}), 'Keys', {'eegid','eegid'});
    groupid_table.subgroups = categorical(strcat(groupid_table.group, '_', groupid_table.sex));
    groupid_table.subidx = grp2idx(groupid_table.subgroups);
    
    grpClinNames = groupid_table.subgroups;
    grpIdxClin = groupid_table.subgroups;
    
clintable_chirp =  {clindata.eegid clindata.group clindata.sex}  
    
%% source statistics

% create dummy sfiles
cfg.predefinedBands = {...
    'delta', '2.5, 4', 'mean'; ...
    'theta', '4.5, 7.5', 'mean';....
    'alpha1', '8, 12', 'mean'; ...
    'alpha2', '10, 12.5', 'mean'; ...
    'beta', '15, 29', 'mean'; ...
    'gamma1', '30, 55', 'mean'; ...
    'gamma2', '65, 90', 'mean'};
cfg.timewindow = [0 80];
cfg.win_length = 2;
cfg.win_overlap = 50;
% create power
% Gather source models
    sFilesMne = sFilesSource2
    
    % Calculate Power per Vertex
    % Dataset B1: Source Absolute Band Power
    sPow.mneAbspowTmp = fx_BstElecPow(sFilesMne, cfg, 'source_abspow');
    %sPow.mneAbspow = fx_BstElecPow(sPow.mneAbspowTmp, cfg, 'source_smooth');
    sPow.mneAbspow = fx_bstReplaceTag(sPow.mneAbspowTmp, 'TEMPLATE');
    
    for i = 1 : length(sPow.mneAbspow)
        
        currentTemplate = sPow.mneAbspow(i);
        curMat = matfile(file_fullpath(currentTemplate.FileName), 'Writable',true);
        curMat.TF = curMat.TF(:,1,1);

        curItpc = itpc(:,i);
        curMat.TF = reshape(curItpc, 15002, 1, 1);
        
        curMat.Freqs = curMat.Freqs(1,:);
        curMat.Freqs = {'chirp', '37,43', 'itpc'};

    end
    

        
cfg.tag              = '';  % defined in each segment 
cfg.timewindow       = [0 80];
cfg.isabs            = 0;
cfg.avgtime          = 1; 
cfg.avgfreq          = 1;
cfg.randomizations   = 2000;
cfg.statisticstype   = 1;
cfg.tail             = 'two';
cfg.correctiontype   = 2;  % 1 default 2 cluster  % FDR manually performed
cfg.minnbchan        = 1;
cfg.clusteralpha     = 0.05;

   % Comparision: FXS vs. TDC, Males Only
cfg.tagstr = 'sFiles_ITPC_All';
sFiles_itc = ...
    process_ft_sourcestatistics_wrapper(...
    sPow.mneAbspow(strcmp(groupids, 'FXS')), ...
    sPow.mneAbspow(strcmp(groupids, 'TDC')), ...
    cfg);

%% Import back into Brainstorm

% Export averaged time and subject cortex view from Brainstorm
% variable name is originalMap

% check ImageGridAmp
originalMap.ImageGridAmp

fxsidx = strcmp('FXS', groupids);
tdcidx = strcmp('TDC', groupids);

itpcres.fxs = mean(itpc(:, fxsidx),2);
itpcres.tdc = mean(itpc(:, tdcidx),2);

newMap = originalMap;
newMap.Comment = 'Mean ITPC 40 Hz FXS';
newMap.ImageGridAmp(:,1) = itpcres.fxs;
newMap.ImageGridAmp(:,2) = itpcres.tdc;

newMap2 = originalMap;
newMap2.Comment = 'Mean ITPC 40 Hz TDC';
newMap2.ImageGridAmp(:,1) = itpcres.tdc;
newMap2.ImageGridAmp(:,2) = itpcres.tdc;

diffres = itpcres.fxs(:,1) - itpcres.tdc(:,1);
newMap3 = originalMap;
newMap3.Comment = 'Diff ITPC 40 Hz FXS-TDC';
newMap3.ImageGridAmp(:,1) = diffres;
newMap3.ImageGridAmp(:,2) = diffres;





