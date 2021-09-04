%% Spectral Events in Murine Data
%  Emperical
%  Author: E. Pedapati (Cincinnati Children's Hospital)
%  Date: 5/11/2021
%  Contact: ernest.pedapati@cchmc.org
%  Note: Enable code folding for sections in MATLAB preferences for
%       readability.

% Outline:
% 01 System: Configuration
% 02a Dataset: Import & Parameters
% 02b Dataset: Transform
% 03a Results: Parameters
% 03b Results: Preallocation
% 03c Results: Computation
% 03d Results: Transform
% 04  Export: Save to CSV/R
%
% Outputs:


%% 01 System: Configuration
% prereqs
% start EEGLAB to setup correct paths to common functions
addpath('C:\Users\ernie\Dropbox\htpminimal');
addpath('C:\Users\ernie\Dropbox\spectralevents');


% system paths
makePath = @(basedir,subfolder) fullfile(basedir, subfolder);
syspath.base        = 'E:/data/se_mea';
syspath.base_cloud  = 'C:/Users/ernie/Dropbox/Papers 2021/GedBounds';
syspath.dataset     = makePath(syspath.base, 'A00_ANALYSIS'); % cleaned data
syspath.script      = makePath(syspath.base_cloud, 'MATLAB');
syspath.class       = makePath(syspath.base_cloud, 'MATLAB/classes');
syspath.external    = makePath(syspath.base_cloud, 'MATLAB/external');
syspath.perms       = makePath(syspath.base, 'perm');
syspath.figures     = makePath(syspath.base, 'figs');
syspath.mat         = makePath(syspath.base, 'mat');
syspath.R           = makePath(syspath.base_cloud, 'R/rawdata');

keyfiles.datacsv    = makePath(syspath.dataset, 'A2105111244_subjTable_MEAStage4.csv');
keyfiles.datamat    = makePath(syspath.dataset, 'A2105111244_subjTable_MEAStage4.mat');

outfile.csv.prep    = makePath(syspath.R, '00_gedB_QI_preprocessing.csv');
outfile.csv.trials  = makePath(syspath.R, '01_gedB_QI_validtrials.csv');

% add paths
cellfun(@(x) addpath(x), {syspath.script; syspath.class; syspath.external} ,'uni', 0);

% end
%% 02a Dataset: Import
%  Parameters
minTrials = 80; % minimum valid trials to include (60 x 2 s = 120 s)

p2 = htpPortableClass;   % MATLAB object / methods / properties
p2.importDataObjects( ...
    keyfiles.datacsv, ...
    keyfiles.datamat, ...
    syspath.dataset );
p2.updateBasePaths( syspath.dataset );



%% Spectral Event Toolbox
p.am = p2;
p.am.setBasePath('E:\data\se_mea');
p.am.updateBasePaths('E:\data\se_mea');
% load trial level data for a single channel
% {p.am.sub(1).EEG.chanlocs.labels}
%%
p.am.htpcfg.logger = log4m.getLogger(fullfile(p.am.htpcfg.scriptPath, 'local/logfiles/', p.am.htpcfg.logfile)); 
p.am.assignConfig2Subjects;

p.am.sub(1).loadDataset('postcomps');
channames = {p.am.sub(1).EEG.chanlocs.labels};
numSubj = length(p.am.sub);

for ci = 1 : 1 % length(channames)
    
    for si = 1 : length(p.am.sub)
        
        p.am.sub(si).loadDataset('postcomps');
        tmpMat = p.am.sub(si).se_exportTrialsByChanName(channames{ci}, 1000);
        % x{si} = tmpMat(:, 1:25);
        x{si} = tmpMat;
        if strcmp(p.am.sub(si).subj_subfolder,'Group1')
            classLabels{si} = 1;
        else
            classLabels{si} = 1;
        end
        p.am.sub(si).unloadDataset;
        
    end
end

%%
studyLabel  = 'MEA_CSD';
bandName    = {'theta', 'alpha1', 'alpha2', 'alpha', 'beta', 'gamma1', 'gamma2'};
bandRange   = {[3.5 7.5], [7.5 10], [10 12.5], [7.5 12.5], [15 30], [30 57], [63 80]};
p.am.sub(1).loadDataset('postcomps');
channames   = {p.am.sub(1).EEG.chanlocs.labels};
numSubj     = length(p.am.sub);
samplesPerTrial = p.am.sub(1).EEG.pnts;
classLabels = {};
maxTrials = 80;
chanDataX = {};

% final data structure
% x  time x trial
qicsv = {};

for si = 1 : length(p.am.sub)
 
    p.am.sub(si).loadDataset('postcomps');
    p.am.sub(si).laplacian_perrinX;
    EEG = p.am.sub(si).EEG;
    p.am.sub(si).storeDataset(EEG, p.am.sub(si).pathdb.signal, p.am.sub(si).subj_subfolder, p.am.sub(si).filename.postcomps);
    for ci = 1 : length(channames)
        
        tmpMat = p.am.sub(si).se_exportTrialsByChanName(channames{ci}, samplesPerTrial);
        tmpMat = tmpMat(:,1:maxTrials);
        chanDataX{ci, si} = tmpMat;
        
    end
    
    qicsv{si, 1} = p.am.sub(si).subj_basename;
    qicsv{si, 2} = p.am.sub(si).subj_subfolder;
    qicsv{si, 3} = p.am.sub(si).EEG.trials;
    qicsv{si, 4} = size(chanDataX,2);
    
    classLabels{si} = qicsv{si, 4};
    p.am.sub(si).unloadDataset;
    
end

% chanDataX 30 channels x 23 subjects
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
        
    parfor ci = 1 : length(channames)
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

%%

for bandi = 1 : length(bandSpectralEvents)
    
    specEv_struct = bandSpectralEvents{bandi};
    
    numSubj = length(specEv_struct); %Number of subjects/sessions
    csvout = {};
    % Event feature probability histograms (see Figure 5 in Shin et al. eLife 2017)
    features = {'eventnumber','maximapowerFOM','duration','Fspan'}; %Fields within specEv_struct
    feature_names = {'event number','event power (FOM)','event duration (ms)','event F-span (Hz)'}; %Full names describing each field
    figure
    colNames = {'id','eegid','group','chan',...
        'eventno_total', 'eventno_mean',...
        'duration_total', 'duration_mean', ...
        'maximapowerFOM_total', 'maximapowerFOM_mean', ...
        'Fspan_total', 'Fspan_mean'};
    for ci = 1 : numel(channames)
        
    for feat_i=1:numel(features)
        feature_agg = [];
        for subj_i=1:numSubj
            
            csvout{subj_i, 1} = bandName{bandi};
            csvout{subj_i, 2} = p.am.sub(subj_i).subj_basename;
            csvout{subj_i, 3} = p.am.sub(subj_i).subj_subfolder;
            csviut(subj_i, 4) = channames{ci};
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
    end
    
    resTable = cell2table(csvout, 'VariableNames', colNames);
    csvFilename = sprintf('se_%s_%s.csv', studyLabel, bandName{bandi});
    outputpath = 'C:\Users\ernie\Dropbox\Papers 2021\GedBounds\R\se\';
    savefile = fullfile(outputpath, csvFilename);
    writetable(resTable, savefile );
    disp(savefile);
    
end


%%
% Set dataset and analysis parameters
eventBand = [30,55]; %Frequency range of spectral events
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
csvFilename = sprintf('se_%s_%s.csv', studyLabel, bandName{1});
outputpath = 'C:\Users\ernie\Dropbox\Papers 2021\GedBounds\R\se\';
savefile = fullfile(outputpath, csvFilename);
writetable(resTable, savefile );