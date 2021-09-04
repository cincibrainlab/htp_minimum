% bst statistics to CSV
% Function library 
% === 
%   - fx_getBstGroupStats         : Get all group statistics results
%   - fx_customSaveResultsCSV     : Save cell array to CSV table
%   - fx_getBstTimeFreqFiles      : Get all subject time frequency results

% post comps to BST (incomplete)

%% identify current study
sStudy = bst_get('Study');
sProtocol = bst_get('ProtocolInfo');
sSubjects =  bst_get('ProtocolSubjects');
sStudyList = bst_get('ProtocolStudies'); 
protocol_name = 'P1_FXS70_TDC71_NOBZO';
atlas = fx_BstGetDKAtlasFromSurfaceMat;
sCortex = in_tess_bst('@default_subject/tess_cortex_pial_low.mat');

pathdb.R = 'C:\Users\ernie\Dropbox\NN_RestingSource\Results\data-raw-revision';
%pathdb.R = '/Volumes/extm1/Dropbox-Maestral/NN_RestingSource/Results/data-raw-revision';

%% Calculate EEG power (Electrode and Source)
sFilesRecordings = bst_process('CallProcess', ...
    'process_select_files_data', [], []);

% check path to raw data, fix if needed.
fx_bstViewInvalidPathToRawRecordings( sFilesRecordings )

% Fix Path if Needed
% correctPath.macmini = '/Volumes/data/P1_73FXS_71_TDC/source';
% correctPath.cometlake = 'Y:\P1_73FXS_71_TDC\source\';
% correctPath.original = '@default_subject/tess_cortex_pial_low.mat';
% oldpath = 'D:\data\P1_73FXS_71_TDC\source\';
% fx_bstFixInvalidPathToRawRecordings( sFilesRecordings, correctPath.cometlake );

%% Process: Power spectrum density (Welch)

% Parameters
cfg.predefinedBands = {'delta', '2.5, 4', 'mean'; ...
    'theta', '4.5, 7.5', 'mean';....
    'alpha1', '8, 12', 'mean'; ...
    'alpha2', '10, 12.5', 'mean'; ...
    'beta', '15, 29', 'mean'; ...
    'gamma1', '30, 55', 'mean'; ...
    'gamma2', '65, 80', 'mean'};
cfg.timewindow = [0 80];
cfg.win_length = 2;
cfg.win_overlap = 50;

sPow.elecAbspow = fx_BstElecPow(sFilesRecordings, cfg, 'abspow');
sPow.elecAbspow = fx_bstAddTag(sPow.elecAbspow, 'Elec_ABSPOW');

sPow.elecAbsNorm = fx_BstElecPow(sFilesRecordings, cfg, 'abspownorm');
sPow.elecAbsNorm = fx_bstAddTag(sPow.elecAbsNorm, 'Elec_ABSNORM');

sPow.elecRelpow = fx_BstElecPow(sFilesRecordings, cfg, 'relpow');
sPow.elecRelpow = fx_bstAddTag(sPow.elecRelpow, 'Elec_RELPOWER');

sPow.elecFFT = fx_BstElecPow(sFilesRecordings, cfg, 'FFT');
sPow.elecFFT = fx_bstAddTag(sPow.elecFFT, 'Elec_FFT');

sPow.elecFFTrel = fx_BstElecPow(sFilesRecordings, cfg, 'FFTrel');
sPow.elecFFTrel = fx_bstAddTag(sPow.elecFFT, 'Elec_FFTrel');

% index of datasets
sPowFields = fieldnames(sPow);

%% extract values by region and write CSV

for sPow_i = 5 : length(sPowFields)
    [sResultsPow.(sPowFields{sPow_i}), VariableNames] = ...
        fx_BstExtractValuesElec(sPow.(sPowFields{sPow_i}));
    outfile_csv = fullfile(pathdb.R, [protocol_name '_' sPowFields{sPow_i} '.csv']);
    fx_customSaveResultsCSV(outfile_csv, sResultsPow.(sPowFields{sPow_i}), VariableNames);
end

sum(in_bst(sPow.elecFFTrel(1).FileName).TF(1,:))
test = sPow.(sPowFields{sPow_i});
in_bst(test(1).FileName).TF(1,1:10)

sResultsPow.(sPowFields{sPow_i})

%%

% get all timefreq files from subjects
sFiles = fx_getBstTimeFreqFiles();

% SUBJECT LEVEL RESULTS
% create subsets depending on analysis get tag from bst)
expected_total_subjects = 141;

    % 1. Relative Power
    sData.relpow = fx_getBstSelectByTag(sFiles, ...
        '| relative | ssmooth3 |', 'select', expected_total_subjects);

    % 2. Absolute Power
    sData.abspow = fx_getBstSelectByTag(sFiles, ...
        'PSD: 79/2000ms Power,FreqBands,ABS | ssmooth2', 'select', expected_total_subjects);

    % 3. Absolute Power: 1/F normalized
    sData.abspownorm = fx_getBstSelectByTag(sFiles, ...
        'PSD: 79/2000ms Power,FreqBands,ABS | ssmooth3 | abs7326 | multiply', 'select', expected_total_subjects);
    
    % 4. Relative Power Electrode Level
        sData.relpowelect = fx_getBstSelectByTag(sFiles, ...
        'NOBZO_NOEPSILON_ELEC', 'select', expected_total_subjects);
    
% index of datasets
sDataFields = fieldnames(sData);

% extract values by region and write CSV
VariableNames = {'id','group',...
    'eegid','label','region','bandname','value'};
for sData_i = 1 : length(sDataFields)
    sResults.(sDataFields{sData_i}) = ...
        fx_BstExtractValues(sData.(sDataFields{sData_i}));
    outfile_csv = fullfile(pathdb.R, [protocol_name '_' sDataFields{sData_i} '.csv']);
    fx_customSaveResultsCSV(outfile_csv, sResults.(sDataFields{sData_i}), VariableNames);
end

%% Brainstorm Source Statistics

% retreive stat result structure
sStats = fx_BstRetrieveGroupStats;
sStatsFields = {sStats.Comment};

% create table of sig. thresholds for each comparision
sStatsSigThreshold = fx_getBstSigThresholds(sStats);
VariableNames = {'Comparison', 'Correction', 'alpha', 'pthreshold'};
outfile_csv = fullfile(pathdb.R, [protocol_name '_sStatsSigThreshold.csv']);
fx_customSaveResultsCSV(outfile_csv, sStatsSigThreshold, VariableNames);

%% Plotting

% get Brainstorm Global Variables
global GlobalData;

availableStats = sStudy.Stat;
orientation = 'fig';
OutputFile = 'C:\Users\ernie\Dropbox\Papers 2021\GedBounds\R\figs';
outputdir = 'C:\Users\ernie\Dropbox\NN_RestingSource\Results\figs_revision';

viewports = {'bottom', 'top',  'right_intern', 'left_intern',  'left', 'right'};

for stati = 1 : length(availableStats)
    OverlayFile = availableStats(stati).FileName;
    OverlayDetails = in_bst(OverlayFile);
    
    
    if ~isempty(OverlayDetails.SurfaceFile)
        [hFig, iDS, iFig] = view_surface_data([], OverlayFile);
        for viewi = 1 : length(viewports)
            currentview = viewports{viewi};
            panel_surface('SelectHemispheres', 'reset');
            figure_3d ('SetStandardView', hFig, currentview)
            for freqi = 1 : length(OverlayDetails.Freqs)
                currentFreq = OverlayDetails.Freqs{freqi, 1};
                GlobalData.UserFrequencies.iCurrentFreq = freqi;
                bst_figures('FireCurrentFreqChanged');
                imgFile = fullfile(outputdir, [OverlayDetails.Comment '_' currentview '_' currentFreq '.tif']);
                
                %[ hContactFig ] = view_contactsheet( hFig, 'freq', orientation, OutputFile);
                tmpimg = out_figure_image( hFig, imgFile, []);
                
                % modify to white background image
                
                % can modify image cropping here
                % RECT is a 4-element vector with the form [XMIN YMIN WIDTH HEIGHT];
                
                function img = fx_customImageMakeWhiteBackground( img )
                [L,W,D] = size(A);
                for i = 1:L
                    for j = 1:W
                        for k = 1:D
                            if A(i,j,k) == 0
                                A(i,j,k) = 255;
                            end
                        end
                    end
                end
                figure; imshow(A);
                imgsize_original = size(rawimg{viewi, freqi, stati});
                imgsize_original = imgsize_original(1:2);
                imgsize = [imgsize(2) imgsize(1)]; % Width, Height
                
                win = centerCropWindow2d(imgsize_original, [450 600])
                figure; imshow(imcrop(rawimg{viewi, freqi, stati}, win));
            rawimg{viewi, freqi, stati} 
            end
        end
        close(hFig);
        
        
        for i = 1 : length(OverlayDetails.Freqs)
            newimg{i, stati} = [rawimg{:, i,stati}];
        end
        imwrite(vertcat(newimg{:, stati}), fullfile(outputdir, [OverlayDetails.Comment '_Panel.tif']))
    end
end
%%
sColormap = db_template('Colormap');
sColormap.MaxMode  = 'custom';
sColormap.MinValue = 6;
sColormap.MaxValue = -6;

fx_BstSetColormap(hFig, sColormap)
    


ColormapInfo = getappdata(hFig, 'Colormap');
        ColormapType = ColormapInfo.Type;

                ColormapInfo = getappdata(hFig, 'Colormap');
backupColormapType = ColormapType;
backupsColormap =  sColormap;

    % Update colormap
    sColormap.MaxMode  = 'custom';
    sColormap.MinValue = 3;
    sColormap.MaxValue = -3;
    SetColormap('stat2', sColormap);
    % Update all figures
    
    
    FireColormapChanged(ColormapType);
    
    
function SetMaxMode(ColormapType, maxmode, DisplayUnits)
    % Parse inputs
    if (nargin < 3) || isempty(DisplayUnits)
        DisplayUnits = [];
    end
    % Check values
    if ~ismember(lower(maxmode), {'local','global','custom'})
        error(['Invalid maximum mode: "' maxmode '"']);
    end
    % Custom: ask for custom values
    if strcmpi(maxmode, 'custom')
        SetMaxCustom(ColormapType, DisplayUnits);
    else
        % Update colormap
        sColormap = GetColormap(ColormapType);
        sColormap.MaxMode = lower(maxmode);
        SetColormap(ColormapType, sColormap);
        % Fire change notificiation to all figures (3DViz and Topography)
        FireColormapChanged(ColormapType);
    end
end


 fullfile(outputdir, [OverlayDetails.Comment '_Panel_' currentFreq '.tif'])
newimage = {};
for celli = 1 : size(rawimg,1)
   
    if celli == 1
        newimag = rawimg{celli,1};
    else
        newimag = [newimage rawimg{celli,1}];
    end
    
end
[rawimg{:,1}]
figure; imshow([rawimg{1,1} rawimg{2,1}]);
figure; imshow([rawimg{:, :,stati}]);
figure; imshow(vertcat(newimg{:, stati}));

stripimg = [rawimg{1} rawimg{3}]
figure; imshow(stripimg)
            figure_3d ('SetStandardView', hFig, 'top')
               panel_surface('SelectHemispheres', 'reset');

% figures: bottom, top, right/split, left/split, right, left
% {'left', 'right', 'top', 'left_intern', 'right_intern', 'bottom'}
% {'bottom', 'top', 'right_intern', 'left_intern', 'right', 'left'}
%%
sStatsData = in_bst(bst_get('AnyFile',sFiles(1).FileName).Stat(1).FileName);
[hFig, iDS, iFig] = view_surface(in_bst(bst_get('AnyFile',sFiles(1).FileName).Stat(1).FileName).SurfaceFile)


            TfInfo = getappdata(hFig, 'Timefreq');
            TfInfo.iFreqs = GlobalData.UserFrequencies.iCurrentFreq;
            setappdata(hFig, 'Timefreq', TfInfo);
            UpdateSurfaceData(hFig);
            



OverlayFile = 'Group_analysis/@intra/ptimefreq_psd_no_200905_1641.mat';
[hFig, iDS, iFig] = view_surface_data(SurfaceFile, OverlayFile, Modality, FigureOption, ShowScouts)
[hFig, iDS, iFig] = view_surface_data([], OverlayFile)

[sStudy, iStudy, iFile, DataType, sFileMat] = bst_get('AnyFile', OverlayFile);
    hFig = GlobalData.DataSet(iDS).Figure(iFig).hFigure;

 TfInfo = getappdata(hFig, 'Timefreq');
 % If there are some time-frequency recordings in this
    if ~isempty(TfInfo)
        % Update frequency to display
        TfInfo.iFreqs = GlobalData.UserFrequencies.iCurrentFreq;
        setappdata(hFig, 'Timefreq', TfInfo);
        % Update display
        UpdateSurfaceData(hFig);
    end

figure_3d ('SetStandardView', hFig, 'top')


%% Report Significant Regions
% get labels counts associated with significant regions, two tails
% including M/F comparisons
resultTable = fx_customStatReportSigLabels( sStats, atlas );

VariableNames = {'statCompare','bandname', 'label','tail1','tail2'};
outfile_csv = fullfile(pathdb.R, [protocol_name '_StatCompareRegions.csv']);
fx_customSaveResultsCSV(outfile_csv, resultTable, VariableNames);

% Report Significant Values per Subject
VariableNames = {'statCompare','eegid', 'group','bandname','tail1', 'tail2'};

resultTable = fx_customStatReportSigPowerPerSubject( sStats(contains(sStatsFields', 'REL_ALL')), sData.relpow );
outfile_csv = fullfile(pathdb.R, [protocol_name '_SigOnly_REL_ALL.csv']);
fx_customSaveResultsCSV(outfile_csv, resultTable, VariableNames);

resultTable = fx_customStatReportSigPowerPerSubject( sStats(contains(sStatsFields', 'REL_Electrode')), sData.relpowelect);
outfile_csv = fullfile(pathdb.R, [protocol_name '_SigOnly_REL_Electrode_ALL.csv']);
fx_customSaveResultsCSV(outfile_csv, resultTable, VariableNames);

resultTable = fx_customStatReportSigPowerPerSubject( sStats(contains(sStatsFields', 'ABS_Raw')), sData.abspow);
outfile_csv = fullfile(pathdb.R, [protocol_name '_SigOnly_ABS_Raw.csv']);
fx_customSaveResultsCSV(outfile_csv, resultTable, VariableNames);

resultTable = fx_customStatReportSigPowerPerSubject( sStats(contains(sStatsFields', 'ABS_Norm')), sData.abspownorm );
outfile_csv = fullfile(pathdb.R, [protocol_name '_SigOnly_ABS_Norm.csv']);
fx_customSaveResultsCSV(outfile_csv, resultTable, VariableNames);

%% get timeseries for advanced analyses

% get Recording Files
sFilesRecordings = bst_process('CallProcess', ...
    'process_select_files_data', [], []);
correctPath.macmini = '/Volumes/data/P1_73FXS_71_TDC/source';
correctPath.cometlake = 'Y:\P1_73FXS_71_TDC\source\';
correctPath.original = '@default_subject/tess_cortex_pial_low.mat';
oldpath = 'D:\data\P1_73FXS_71_TDC\source\';

fx_bstFixInvalidPathToRawRecordings( sFilesRecordings, correctPath.cometlake );

% can create source or electrode time series at this juncture


sFilesSource = fx_getBstSourceFiles();

for sourcei = 1 : length(sFilesSource)
sFilesTS = bst_process('CallProcess', 'process_extract_scout', sFilesSource(sourcei), [], ...
                        'scouts',         {atlas.Name, {atlas.Scouts.Label}}, ...
                        'scoutfunc',      1, ...  % PCA
                        'isflip',         1, ...
                        'isnorm',         0, ...
                        'concatenate',    1, ...
                        'save',           1, ...
                        'addrowcomment',  1, ...
                        'addfilecomment', 1);
end                

sFilesMatrix = fx_getBstMatrixFiles();
 
sData.TS = fx_getBstSelectByTag(sFilesMatrix, ...
        'Tag_82024', 'select', expected_total_subjects);
    
sMat = in_bst(sData.TS(1).FileName); % contains raw values

% sProps = in_bst_data(sData.TS(1).FileName); % contains properies, srate etc. 0.1158 -0.5111

sMat2T = in_bst(sFilesRecordings(1).FileName);
sMat2P = in_bst_data(sFilesRecordings(1).FileName);
bst_get('ChannelFile', '@default_study/channel.mat')

  function corr_rho = SpearmanCorrelation2(rel_target_band, rel_gamma_band)
            [rho, ~] = corr(rel_target_band, rel_gamma_band, 'Type', 'Spearman');
            corr_rho = mean(rho);
        end
        
        
        function corr_rho = SpearmanCorrelation(rel_target_band, rel_gamma_band)
            cluster_value = mean(rel_target_band, 2); % 1* trials
            [rho, ~] = corr(repmat(cluster_value, [1 size(rel_gamma_band,2)]), ...
                                rel_gamma_band, 'Type', 'Spearman');
            corr_rho = rho(1, :);
        end

        
        %% create Grp Ids with Sex
        
[eegids, groupids] = fx_customGetSubNames(sData.relpow);
        
% add M/F 
csvfile = 'C:\Users\ernie\Dropbox\NN_RestingSource\Results\data-raw\fxssource_table_01_demographics_data.csv';
[subgroups subgroup_idx] = fx_customAddClinicalVarByEegid( sData.relpow, csvfile);
        