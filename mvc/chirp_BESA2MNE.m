% chirp BESA generic -> Brainstorm MNE -> EEG SET
% convert to EEG SET for import to Brainstorm
% 4/25/2021

%% common import files
cfg = struct();
cfg.chanlocs_file = 'chanfiles/GSN-HydroCel-129_new.sfp';
cfg.import_directory = 'E:/onedrive/OneDrive - cchmc/Datashare - EP LAB/EEG/Chirp_ST/';
cfg.export_directory = 'E:/data/mne_chirp/S04_POSTCOMPS';
cfg.basePath = 'E:/data/mne_chirp/';
cfg.condition = 'Chirp';
cfg.comments = 'EEGSET Pre-MNE';
cfg.study_title = 'MNECHIRP';
cfg.points_per_trial = 1626;
cfg.srate = 500;
cfg.xmin  = -.5;
cfg.number_channels = 129;

am = htpPreprocessMaster;
am.firstRun;
am.preprocess_direct(cfg);

%%
% perform source modeling
cfg.basePath = 'E:/data/mne_chirp/';
datapath = cfg.basePath;
csvfile_subjects = 'A2104252057_subjTable_htp2bst_Stage4.csv';
chanlocsMat = '128_chanlocs.mat';

p = publishClass;
pcfg = [];
p.assignDataPath( datapath );
p.assignSubjectCsv( csvfile_subjects );
p.assignChanLocsMat( chanlocsMat );

%p.assignInitialParameters( pcfg );
p.createAnalysisMaster( true );  % use 'true' to refresh object
p.am.htpcfg.chanNow = p.am.htpcfg.chaninfo(9);

p.am.updateBasePaths( p.pcfg.datapath );
p.am.sub(1).createPaths;
grpInfo = p.am.getGroupInfo;
% p.exportStudyVariables;
p.am.assignChanInfoToSubjects( p.am.getSubjectNetName );
p.am.assignConfig2Subjects;
%%
brainstorm;

p.am.bst_createNewProtocol(...
    'MNEChirp8');
p.am.bst_prepareContinuousData;
p.am.bst_controller_htp('addSubjectToProtocol_rest');

chanInfoStruct.headModel = 'ICBM152';
chanInfoStruct.brand = 'GSN';
chanInfoStruct.chanNumber = '128';
chanInfoStruct.chanLabelFormat = 'E1';

[netIndex, allNetOptions] = p.am.bst_locateChannels(chanInfoStruct);

p.am.bst_controller_recordings('fixChannels');
p.am.bst_reloadCurrentProtocol;
p.am.bst_controller_recordings('computeHeadModel');
p.am.bst_controller_recordings('computeIdentityNoisecov');
p.am.bst_reloadCurrentProtocol;

p.am.bst_computeSources('mne_default');

% open surface file in GUI
sFiles = p.am.bst_getAllSources;

%% export TS
cfg = [];
%cfg.timew = [0 80];
cfg.atlasname = 'Mindboggle';
cfg.atlasname = 'Desikan-Killiany';
cfg.saveSignalSets = true;
cfg.calcPowerTable = false; % suppress calPowerTable in 'generate_PowTable'

%%% check-in BrainStorm
atlas = p.am.bst_getAtlas;
% Name: 'Desikan-Killiany'
% Scouts: [1×68 struct]
atlas = atlas(3);
cfg.atlasname = 'Desikan-Killiany';
cfg.scouts = {atlas.Scouts.Label};
cfg.atlas = atlas;
res = p.am.extractScoutsTimeSeriesPCA(cfg, sFiles,  grpInfo );
%res = p.am.generate_PowTable(cfg, sFiles,  grpInfo ); % export signal ONLY
p.am.createResultsCsv( p.am.sub, 'signal', 'signal' )

% export dataset times series


%% new CSV loop
studyLabel  = 'CHIRP';
bandName    = {'theta', 'alpha1', 'alpha2', 'alpha', 'beta', 'gamma1', 'gamma2'};
bandRange   = {[3.5 7.5], [7.5 10], [10 12.5], [7.5 12.5], [15 30], [30 57], [63 80]};
p.am.sub(1).loadDataset('signal');
channames   = {p.am.sub(1).EEG.chanlocs.labels};
channames   = {'superiortemporalL','superiortemporalR', 'banksstsL','banksstsR', 'transversetemporalR',  'transversetemporalL', ...
    'cuneusL', 'cuneusR'};
numSubj     = length(p.am.sub);
samplesPerTrial = 1626;
classLabels = {};
maxTrials = 80;
chanDataX = {};
%% generate publication figure
searchSub_FXS = 'D0320_chirp-ST';
searchSub_TDC = 'D1486_chirp-ST';
searchsub = searchSub_FXS;
find(strcmp(searchsub, {p.am.sub.subj_basename}))

studyLabel  = 'CHIRP';
bandName    = {'gamma1' };
bandRange   = {[30 55]};
p.am.sub(1).loadDataset('signal');
channames   = {p.am.sub(1).EEG.chanlocs.labels};
channames   = {'superiortemporalR'};
numSubj     = length(p.am.sub);
samplesPerTrial = 1626;
classLabels = {};
maxTrials = 80;
chanDataX = {};

%%
%% generate MEA publication figure
%searchSub_FXS = 'D0320_chirp-ST';
%searchSub_TDC = 'D1486_chirp-ST';
%searchsub = searchSub_TDC;
%find(strcmp(searchsub, {p.am.sub.subj_basename}))

studyLabel  = 'MEA_REST';
bandName    = {'gamma1' };
bandRange   = {[30 55]};
p.am.sub(1).loadDataset('signal');
channames   = {p.am.sub(1).EEG.chanlocs.labels};
channames   = {'Ch 12'};
numSubj     = length(p.am.sub);
samplesPerTrial = 1000;
classLabels = {};
maxTrials = 80;
chanDataX = {};


%%
% final data structure
% x  time x trial


%%
qicsv = {};
si_count = 1;
classLabels = {};
chanDataX = {};
for si = 14:14 %1: length(p.am.sub)

    s = p.am.sub(si);
    s.loadDataset('signal');
    % p.am.sub(si).laplacian_perrinX;
    for ci = 1 : length(channames)
        
        tmpMat = s.se_exportTrialsByChanName(channames{ci}, samplesPerTrial);
       % tmpMat = tmpMat(:,1:maxTrials);
        chanDataX{ci, si_count} = tmpMat;
        
    end
    
    qicsv{si_count, 1} = s.subj_basename;
    qicsv{si_count, 2} = s.subj_subfolder;
    qicsv{si_count, 3} = s.EEG.trials;
    qicsv{si_count, 4} = size(chanDataX{ci, si_count},2);
    
    classLabels{si_count} = qicsv{si_count, 4};
    s.unloadDataset;
    
    si_count = si_count + 1;
end

%%
B =reshape(x{1}, 1,[]);
C = reshape(B(1:1626*81), 1626,[]);
x={C};
%%
    chanSpectralEvents = cell(numel(channames),1);
    chanTFRs = cell(numel(channames),1);
    chantimeseries = cell(numel(channames),1);
    
bandSpectralEvents = {};
bandSpectralEvents = {};
bandTFRs = {};
bandtimeseries = {};


for bi = 1 : length(bandName)
        disp(bandName{bi});
        disp(bandRange{bi});
        
    for ci = 1 : length(channames)
        disp(['Channel:' channames{ci}]);     
        
        eventBand = bandRange{bi};
        fVec = (2:80);
        Fs = 500;
        
        findMethod = 1; %Event-finding method (1 allows for maximal overlap while 2 limits overlap in each respective suprathreshold region)
        vis = true;
        
        x = chanDataX(ci, :);
        
        [specEvents,TFRs,timeseries] = ...
            spectralevents(eventBand,fVec, ...
            Fs,findMethod,vis,x, ...
            classLabels); %Run spectral event analysis
       
    chanSpectralEvents{ci} = specEvents;
    chanTFRs{ci} = TFRs;
    chantimeseries{ci} = timeseries;  
    
    end
    
    % bandSpectralEvents{bi} = chanSpectralEvents;
    % bandTFRs{bi} = chanTFRs;
    %bandtimeseries{bi} = chantimeseries;
    csvout = [];
    tcount = 1;
    colNames = {'id','eegid','group','chan',...
        'eventno_total', 'eventno_mean',...
        'duration_total', 'duration_mean', ...
        'maximapowerFOM_total', 'maximapowerFOM_mean', ...
        'Fspan_total', 'Fspan_mean'};
    
    for ci = 1 : length(chanSpectralEvents)
        
        specEvents = chanSpectralEvents{ci};
        specEv_struct = specEvents;
        
        numSubj = length(specEv_struct); %Number of subjects/sessions
        
        % Event feature probability histograms (see Figure 5 in Shin et al. eLife 2017)
        features = {'eventnumber','maximapowerFOM','duration','Fspan'}; %Fields within specEv_struct
        feature_names = {'event number','event power (FOM)','event duration (ms)','event F-span (Hz)'}; %Full names describing each field
        
        for feat_i=1:numel(features)
            feature_agg = [];
            for subj_i=1:numSubj
                
                csvout{subj_i, 1} = bandName{bi};
                csvout{subj_i, 2} = p.am.sub(subj_i).subj_basename;
                csvout{subj_i, 3} = p.am.sub(subj_i).subj_subfolder;
                csvout{subj_i, 4} = channames{ci};
                % Feature-specific considerations
                if isequal(features{feat_i},'eventnumber')
                    feature_agg = [feature_agg; ...
                        specEv_struct(subj_i).TrialSummary.TrialSummary.(features{feat_i})];
                    csvout{subj_i, 5} = sum(specEv_struct(subj_i).TrialSummary.TrialSummary.(features{feat_i}));
                    csvout{subj_i, 6} = mean(specEv_struct(subj_i).TrialSummary.TrialSummary.(features{feat_i}));
                    
                else
                    if isequal(features{feat_i},'duration')
                        feature_agg = [feature_agg; specEv_struct(subj_i).Events.Events.(features{feat_i}) * 1000]; %Note: convert from s->ms
                        csvout{subj_i, 7} = sum(specEv_struct(subj_i).Events.Events.(features{feat_i}) * 1000);
                        csvout{subj_i, 8} = mean(specEv_struct(subj_i).Events.Events.(features{feat_i}) * 1000);
                        
                    else
                        if isequal(features{feat_i},'maximapowerFOM')
                            feature_agg = [feature_agg; specEv_struct(subj_i).Events.Events.(features{feat_i})];
                            csvout{subj_i, 9} = sum(specEv_struct(subj_i).Events.Events.(features{feat_i}) * 1000);
                            csvout{subj_i, 10} = mean(specEv_struct(subj_i).Events.Events.(features{feat_i}) * 1000);
                            
                        else
                            if  isequal(features{feat_i},'Fspan')
                                feature_agg = [feature_agg; specEv_struct(subj_i).Events.Events.(features{feat_i})];
                                csvout{subj_i, 11} = sum(specEv_struct(subj_i).Events.Events.(features{feat_i}) * 1000);
                                csvout{subj_i, 12} = mean(specEv_struct(subj_i).Events.Events.(features{feat_i}) * 1000);
                                
                            else
                                
                            end
                        end
                    end
                end
                
                
            end
            
        end
        
        if ci == 1
            csvtotal = csvout;
        else
        csvtotal = [csvtotal ; csvout];
        end
    end
    
    resTable = cell2table(csvtotal, 'VariableNames', colNames);
    csvFilename = sprintf('se_%s_%s.csv', studyLabel, bandName{bi});
    outputpath = 'C:\Users\ernie\Dropbox\Papers 2021\GedBounds\R\se\';
    savefile = fullfile(outputpath, csvFilename);
    writetable(resTable, savefile );
    disp(savefile);
    
    
end


%% Spectral Event Toolbox

% load trial level data for a single channel
% {p.am.sub(1).EEG.chanlocs.labels}

for si = 1 : length(p.am.sub)
    
    p.am.sub(si).loadDataset('signal');
    tmpMat = p.am.sub(si).se_exportTrialsByChanName('superiortemporalR', 1626);
    % x{si} = tmpMat(:, 1:25);
    x{si} = tmpMat;
    if strcmp(p.am.sub(si).subj_subfolder,'FXS')
        classLabels{si} = 1;
    else
        classLabels{si} = 1;
    end
    p.am.sub(si).unloadDataset;
    
end


%%
% Load data sessions/subjects from the same experimental setup so that 
% spectral event features are differentially characterized only between the 
% desired trial classification labels: in this case, detection vs. 
% non-detection
numSubj = 10;
x = cell(1,numSubj);
classLabels = cell(1,numSubj);
for subj_i=1:numSubj
    load(['test_data/prestim_humandetection_600hzMEG_subject',num2str(subj_i),'.mat'])
    x{subj_i} = prestim_raw_yes_no'; %Time-by-trial matrix of timeseries trials for detection/non-detection prestimulus MEG
    classLabels{subj_i} = YorN'; %Column vector of trial classification labels
    clear prestim_raw_yes_no YorN
end
%%
% Set dataset and analysis parameters
eventBand = [63,80]; %Frequency range of spectral events
fVec = (2:80); %Vector of fequency values over which to calculate TFR
Fs = 500; %Sampling rate of time-series
findMethod = 1; %Event-finding method (1 allows for maximal overlap while 2 limits overlap in each respective suprathreshold region)
vis = false; %Generate standard visualization plots for event features across all subjects/sessions
%tVec = (1/Fs:1/Fs:1);
[specEvents,TFRs,timeseries] = ...
    spectralevents(eventBand,fVec,Fs,findMethod,vis,x,classLabels); %Run spectral event analysis

%%
spectralevents_vis2(specEvents)


% Save figures
classes = [0,1];
for subj_i=1:numSubj
    for class_i=1:2
        figName = strcat('./test_results/matlab/prestim_humandetection_600hzMEG_subject', num2str(subj_i), '_class_', num2str(classes(class_i)), '.png');
        saveas(figure((subj_i-1)*2+class_i),figName);
    end
end

specEvents

%%
specEv_struct = specEvents;

numSubj = length(specEv_struct); %Number of subjects/sessions
csvout = {};
% Event feature probability histograms (see Figure 5 in Shin et al. eLife 2017)
features = {'eventnumber','maximapowerFOM','duration','Fspan'}; %Fields within specEv_struct
feature_names = {'event number','event power (FOM)','event duration (ms)','event F-span (Hz)'}; %Full names describing each field
figure
colNames = {'id','eegid','group',...
    'eventno_total', 'eventno_mean',...
    'duration_total', 'duration_mean', ...
    'maximapowerFOM_total', 'maximapowerFOM_mean', ...
    'Fspan_total', 'Fspan_mean'};
for feat_i=1:numel(features)
    feature_agg = [];
    for subj_i=1:numSubj
        
        csvout{subj_i, 1} = subj_i;
        csvout{subj_i, 2} = p.am.sub(subj_i).subj_basename;
        csvout{subj_i, 3} = p.am.sub(subj_i).subj_subfolder;

        % Feature-specific considerations
        if isequal(features{feat_i},'eventnumber')
            feature_agg = [feature_agg; ...
                specEv_struct(subj_i).TrialSummary.TrialSummary.(features{feat_i})];
            csvout{subj_i, 4} = sum(specEv_struct(subj_i).TrialSummary.TrialSummary.(features{feat_i}));
            csvout{subj_i, 5} = mean(specEv_struct(subj_i).TrialSummary.TrialSummary.(features{feat_i}));
            
        else
            if isequal(features{feat_i},'duration')
                feature_agg = [feature_agg; specEv_struct(subj_i).Events.Events.(features{feat_i}) * 1000]; %Note: convert from s->ms
                csvout{subj_i, 6} = sum(specEv_struct(subj_i).Events.Events.(features{feat_i}) * 1000);
                csvout{subj_i, 7} = mean(specEv_struct(subj_i).Events.Events.(features{feat_i}) * 1000);
                
            else
                if isequal(features{feat_i},'maximapowerFOM')
                    feature_agg = [feature_agg; specEv_struct(subj_i).Events.Events.(features{feat_i})];
                    csvout{subj_i, 8} = sum(specEv_struct(subj_i).Events.Events.(features{feat_i}) * 1000);
                    csvout{subj_i, 9} = mean(specEv_struct(subj_i).Events.Events.(features{feat_i}) * 1000);
                    
                else
                    if  isequal(features{feat_i},'Fspan')
                        feature_agg = [feature_agg; specEv_struct(subj_i).Events.Events.(features{feat_i})];
                        csvout{subj_i, 10} = sum(specEv_struct(subj_i).Events.Events.(features{feat_i}) * 1000);
                        csvout{subj_i, 11} = mean(specEv_struct(subj_i).Events.Events.(features{feat_i}) * 1000);
                        
                    else
                    end
                end
            end
        end
    end
end   

resTable = cell2table(csvout, 'VariableNames', colNames);
writetable(resTable, 'C:\Users\ernie\Dropbox\Papers 2021\GedBounds\R\se\se_higamma.csv');