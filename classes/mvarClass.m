classdef mvarClass < handle
    %mvarClass
    %   Source separation using gneralized eigendecomposition
    %   Code sections adapted from Marrit Zuure & Michael X Cohen
    %   AUTHOR: ernest.pedapati@cchmc.org
    properties
        s   % eegDataClass subject level class
        EEGpresent;
        isSegmentedData;
        lastCovMat;
        lastCovdistz;
        lastEEG;
        
        pathdb; 
        
        permutation_path;
        figs_path;
        
        currentSampleRate;
        
        SCovMat;
        RCovMat;
        SCovComp;
        SCovCompNumber;
        RCovComp;
        CompTimes;
        atlas;
        filtData;
        evecsT;
        evals;
        maxcomp;
        lastFilterHz
        lastFilterFhwm
        lastFilterKernal
        lastFilterOrder
        lastFilterLowerBound
        lastFilterUpperBound
        saveplot
        lastEvecs
        lastEvals
        
        lastGedResults;
        lastMidF;
        
        filedb;
        
        
        covS_window
        covR_window
        
        permresults
        
        % PAC functions
        pac_covT;
        
        % gedbounds
        gedBoundsRes
        
        minTrials
        lowerBoundGED
        upperBoundGED
        freqStepsGED
        
        topoRow
        
    end
    methods
        function obj = mvarClass()
            %MVARCLASS Construct an instance of this class
            %
            obj.EEGpresent = false;
            obj.saveplot = false;
        end
        % load subject
    end
    methods % housekeeping functions
        function obj = loadSubject(obj, s, datatype)
            % most commonly 'signal' (MNE) or 'postcomps' (cleaned)
            if nargin < 3
                dataLoad = false;
            else
                dataLoad = true;
            end
            
            try
                obj.s = s;
                if dataLoad == true
                    obj.s.loadDataset(datatype);
                    obj.s.EEG.data = double(obj.s.EEG.data);
                    obj.currentSampleRate = obj.s.EEG.srate;
                    obj.EEGpresent = true;
                    obj.getEEGstats;
                else
                    obj.EEGpresent = false;
                end
            catch
                disp('mvarClass: Data loading error.');
            end
        end
        function obj = setPathDb( obj, pathdb)
            obj.pathdb = pathdb;
        end
        function pathdb = getPathDb( obj )
            pathdb = obj.pathdb; 
        end
        function obj = manualLoadEegFile(obj, filename)
            obj.s = eegDataClass;
            obj.s.EEG = pop_loadset(filename);
            obj.s.EEG.data = double(obj.s.EEG.data);
            obj.EEGpresent = true;
            obj.currentSampleRate = obj.s.EEG.srate;
            obj.getEEGstats;
        end
        function obj = setCurrentEEG( obj, EEG )
            obj.s.EEG = EEG;
            obj.getEEGstats;
        end
        function obj = setPermutationDirectory( obj, filepath )
            filepath = fullfile(filepath);
            obj.permutation_path = filepath;
        end
        function obj = setFigsDirectory( obj, filepath )
            filepath = fullfile(filepath);
            obj.figs_path = filepath;
        end
        function filepath  = getPermutationDirectory(obj )
            filepath = obj.permutation_path;
        end
        function filepath = getFigsDirectory(obj)
             filepath = obj.figs_path;
        end
        function obj = keepScalpElectrodesOnly( obj )
            EEG = obj.s.EEG;
            keepElectrodes = obj.atlas.labels;
            EEG = pop_select(EEG, ...
                'channel', ...
                keepElectrodes); %scalp electrode
            obj.s.EEG = EEG;
            fprintf('\n\n%d scalp electrodes selected.\n\n', ...
                height(obj.atlas));
        end
        function obj = unloadSubject(obj)
            obj.s = [];
            obj.EEGpresent = false;
            disp('mvarClass: Data unloaded.');
        end
        function obj = getEEGstats( obj )
            if obj.EEGpresent == true
                EEG = obj.s.EEG;
                if EEG.trials > 1
                    trialmsg = sprintf('TRIAL DATA (TRIALS=%d)', EEG.trials);
                    obj.isSegmentedData = true;
                else
                    trialmsg = 'CONTINUOUS DATA (TRIALS=1)';
                    obj.isSegmentedData = false;
                end
                fprintf('\n-----------------------\n');
                fprintf('EEG data (%s)\n', EEG.subject);
                fprintf('-----------------------\n');
                fprintf('%s\n', trialmsg);
                fprintf('-----------------------\n');
                fprintf('Group: %s\n', EEG.group);
                fprintf('Nbchan: %d\t', EEG.nbchan);
                fprintf('Srate: %d\n', EEG.srate);
                fprintf('Data: %s\n', mat2str(size(EEG.data)));
                fprintf('xmin: %1.2f\n', EEG.xmin);
                fprintf('xmax: %1.2f\n', EEG.xmax);
                fprintf('Points: %d\n', EEG.pnts);
                try
                    fprintf('Events: %s\n', strjoin(unique({obj.s.EEG.event.type})));
                catch
                    fprintf('Data not segmented. \n');
                end
                fprintf('-----------------------\n');
                fprintf('to view data use: ged.plotDataTimeSeries\n\n');
            else
                disp('\n\getEEGstats (nmvarClass): No data loaded.\n\n');
            end
            
        end
        function obj = epoch2continuous( obj )
            
            % resting data only as it does not maintain old events
            EEG = obj.s.EEG;
            EEG.data = reshape(EEG.data, size(EEG.data,1), size(EEG.data,2)*size(EEG.data,3));
            EEG.pnts   = size(EEG.data,2);
            EEG.trials = 1;
            EEG.event = [];
            EEG.epoch = [];
            EEG.urevent = [];
            EEG.times  = [];
            
            obj.s.EEG = eeg_checkset( EEG );
            
            obj.getEEGstats;
        end
        function obj = createRegEpochs( obj, arg_recurrence, arg_limits )
            
            obj.s.EEG = ...
                eeg_regepochs( obj.s.EEG, ...
                'recurrence', arg_recurrence,...
                'limits', arg_limits,...
                'extractepochs', 'on',...
                'eventtype', 'CONT', ...
                'rmbase', NaN );
            obj.s.EEG = eeg_checkset(obj.s.EEG);
            obj.getEEGstats;
            
        end
        function obj = listEvents( obj )
            disp(unique({obj.s.EEG.event.type}));
        end
        function obj = createEventEpochs( obj, eventmarker, limits)
            obj.s.EEG = pop_epoch(...
                obj.s.EEG, eventmarker, ...
                limits, ...
                'verbose', 'yes', ...
                'epochinfo', 'yes');
            obj.getEEGstats;
        end
        function obj = clearLastEEG( obj )
            obj.lastEEG = [];
        end
        function obj = loadSensorInfo( obj, type )
            switch type
                case 'egi128'
                    fprintf('\n\nLoading atlas_dk_networks_mni.csv\nAccess with method: getAtlas\n\n');
                    obj.atlas = readtable('egi128wregions.csv');
                case 'dkatlas'
                    fprintf('\n\nLoading atlas_dk_networks_mni.csv\nAccess with method: getAtlas\n\n');
                    obj.atlas = readtable('atlas_dk_networks_mni.csv');
                otherwise
                    fprintf('\n\nloadSensorInfo (mvarClass): Please define %s\n\n', type);
            end
        end
        function atlas = getAtlas( obj )
            atlas = obj.atlas;
        end
        function obj = setSCovMat( obj, covmat )
            obj.SCovMat = covmat;
        end
        function obj = setRCovMat( obj, covmat )
            obj.RCovMat = covmat;
        end
        function EEG = getFiltData( obj )
            EEG = obj.filtData;
        end
        function obj = SavePlotsOn( obj )
            obj.saveplot = true;
        end
        function obj = SavePlotsOff( obj )
            obj.saveplot = false;
        end
        function EEG = getCurrentEEG( obj )
            EEG = obj.s.EEG;
        end
        function eegid = getEegId( obj )
            eegid = obj.s.subj_basename;
        end
        function res = getGroup( obj )
            res = obj.s.subj_subfolder;
        end
        function ged = getLastGedResults( obj )
            ged = obj.lastGedResults;
        end
        function obj = selectEEGTimeRange( obj, num_of_seconds )
            EEG = obj.getCurrentEEG;
            assert(obj.isSegmentedData == false, ...
                'Data is not continuous. Use trial selection instead.');
            EEG = pop_select(EEG, 'time', [1 num_of_seconds] );
            
            obj.setCurrentEEG(EEG);
            
        end
        function obj = reloadCurrentEEG( obj, datatype )
            if nargin < 2
                datatype = 'postcomps';
            end
            obj.s.loadDataset(datatype);
            obj.s.EEG.data = double(obj.s.EEG.data);
            obj.currentSampleRate = obj.s.EEG.srate;
            obj.EEGpresent = true;
            obj.getEEGstats;
        end
        function obj = selectConsecutiveTrials(obj, num_of_trials )
            % num_of_trials= 30;
            EEG = obj.getCurrentEEG;
            
            assert(num_of_trials <= EEG.trials, 'Error: Requested More Trials than Available.')
            selTrialIdx = 1 : num_of_trials;
            EEG = pop_select(EEG, 'trial', selTrialIdx );
            EEG.etc.selTrialGed = selTrialIdx;
            
            obj.setCurrentEEG(EEG);
        end
        function obj = selectRandomTrials( obj, num_of_trials )
            
            % num_of_trials= 30;
            EEG = obj.getCurrentEEG;
            
            assert(num_of_trials <= EEG.trials, 'Error: Requested More Trials than Available.')
            selTrialIdx = randperm(EEG.trials, num_of_trials);
            EEG = pop_select(EEG, 'trial', selTrialIdx );
            EEG.etc.selTrialGed = selTrialIdx;
            
            obj.setCurrentEEG(EEG);
            
        end
        function res = getGedBoundsRes( obj )
            res = obj.gedBoundsRes;
        end
        function res = getTotalTrials( obj )
            EEG = obj.getCurrentEEG;
            res = EEG.trials;
        end
        function res = areRequiredTrialsPresent( obj, desired_trials )
            if obj.getTotalTrials >= desired_trials, res = true; else, res= false; end
        end
        function res = getBadCompNo( obj )
            res = length(obj.s.proc_removeComps);
        end
        function res = getBadChanNo( obj )
            res = length(obj.s.proc_badchans);
        end
        function obj = setMinTrials( obj, minTrials)
            obj.minTrials = minTrials;
        end
        function obj = setGedParameters( obj, lowerBoundGED, upperBoundGED, freqStepsGED)
            obj.lowerBoundGED = lowerBoundGED; 
            obj.upperBoundGED = upperBoundGED;
            obj.freqStepsGED = freqStepsGED;
        end
    end
    methods % computational functions
        function covmat = createCovAvgTrials( obj, inputEEG )
            if nargin < 2
                EEG = obj.s.EEG;
            else
                EEG = inputEEG;
            end
            
            if obj.isSegmentedData == true
                
                % compute average of single-trial covariances
                covave = zeros( EEG.nbchan );
                for triali=1:EEG.trials
                    covave = covave + cov( double(squeeze(EEG.data(:,:,triali))') );
                end
                
                % divide by number of trials
                covave = covave / triali;
                covmat = covave;
                
                obj.lastCovMat = covmat;
                disp('Created Covariance Matrix Averaged from Trials.');
                fprintf('Size (Check!): %s\n', mat2str(size(covmat)));
                fprintf('Min: %e\t', min(max(covmat)));
                fprintf('Max: %e\n', max(max(covmat)));
                
                disp('To visualize use method: plotCovMat');
            else
                fprintf('\n\ncreateCovAvgTrials (mvarClass): No trials found.\n\n');
                obj.getEEGstats;
            end
            
        end
        function covmat = createCovContinous( obj )
            
            if obj.isSegmentedData == false
                EEG = obj.s.EEG;
                % compute average of single-trial covariances
                covave_check = zeros( EEG.nbchan );
                covave = cov(double(obj.s.EEG.data'));
                
                if size(covave_check) ~= size(covave)
                    error('Error: Not Channel Matrix');
                end
                
                obj.lastCovMat = covave;
                covmat = covave;
                disp('Created Covariance Matrix from Continous data.');
                fprintf('Size (Check for channo x channo!): %s\n', mat2str(size(covmat)));
                fprintf('Min: %e\t', min(max(covmat)));
                fprintf('Max: %e\n', max(max(covmat)));
                
                disp('To visualize use method: plotCovMat');
                
            else
                fprintf('\n\ncreateCovAvgTrials (mvarClass): Trial data found.\n\n');
                obj.getEEGstats;
            end
            
        end
        function obj = cleanTrialsByDistToMean(obj, covmat, thresh)
            if nargin < 3
                thresh = 2.3; % ~.01
            end
            
            if obj.isSegmentedData == true
                EEG = obj.s.EEG;
                obj.lastEEG = EEG;
                fprintf('\ncleanTrailsByDistToMean (mvarClass): Last EEG created.\n')
                covave = covmat;
                % now loop through trials and compute the distance to the average
                covdist = zeros(EEG.trials,1);
                
                for triali=1:EEG.trials
                    thistrialcov = cov( squeeze(EEG.data(:,:,triali))' );
                    
                    % compute Frobenius distance
                    covdist(triali) = sqrt( sum(thistrialcov(:) .* covave(:)) );
                    % previous line is the same as this one:
                    %covdist(triali) = sqrt( trace(thistrialcov*covave) );
                    
                    % alternative: Euclidean distance (gives similiar results)
                    %covdist(triali) = sqrt( sum((thistrialcov(:) - covave(:)).^2) );
                end
                
                % convert to z
                covdistz = (covdist-mean(covdist)) / std(covdist);
                obj.lastCovdistz = covdistz;
                
                % threshold
                % defined input
                
                % identify trials that exceed the threshold
                toofar = covdistz>thresh;
                
                % remove those trials from the data
                EEG = pop_select(EEG, 'notrial', toofar);
                fprintf('\nReject trials based on Distance to Mean');
                fprintf('\nTotal rejected Trials: %d.\n\n', sum(toofar));
                obj.s.EEG = EEG;
                obj.getEEGstats;
                
            else
                fprintf('\n\ncreateCovAvgTrials (mvarClass): No trials found.\n\n');
                obj.getEEGstats;
            end
            
        end
        function obj = filterFGx_Cohen( obj, f, fwhm )
            obj.lastFilterHz = f;
            obj.lastFilterFhwm = fwhm;
            EEG = obj.s.EEG;
            obj.filtData = EEG;
            obj.filtData.data = filterFGx(EEG.data,EEG.srate,f,fwhm);
            fprintf('\n\nSignal Filtered @ %d Hz (fwhm: %d)', f, fwhm);
            fprintf('\nTo see signal use method plotLastFiltTimeSeries.\n\n');
        end
        function obj = filtfilt( obj, f, fwhm, EEG )
            if nargin > 3
                % use input EEG
            else
                EEG = obj.s.EEG;
            end
            obj.lastFilterHz = f;
            obj.lastFilterFhwm = fwhm;
            
            %obj.filtData = EEG;
            obj.filtData = ...
                pop_eegfiltnew( EEG, f-fwhm, f+fwhm, [], 0 );
            %eegplot(obj.filtData.data)
            %filterFGx(EEG.data,EEG.srate,f,fwhm);
            fprintf('\n\nSignal Filtered @ %2.2f Hz (fwhm: %2.2f)', f, fwhm);
            fprintf('\nTo see signal use method plotLastFiltTimeSeries.\n\n');
        end
        function obj = perform_ged_wrapper( obj, num_perm )
            % wrapper function for PERFORMGED by  Marrit Zuure
            % available on github
            
            % get filter kernal
            covS_window = obj.covS_window;
            covR_window = obj.covR_window;
            kernel = obj.lastFilterKernal;
            EEG = obj.getCurrentEEG;
            
            [evals, evecs, evecs_rel, covS, covR, evecs_unnormed] = ...
                performGED(EEG, covS_window, covR_window, kernel);
            orig_evals = evals;
            orig_evecs = evecs_unnormed;
            
            obj.perform_gedperm_wrapper( num_perm );
            
            %% Determine significant components
            num_comps = sum(evals > obj.permresults.perm_evals');
            
            %% Cleanup: Drop components that potentially represent eigenplanes
            
            % Flag significant eigenvalues that are closer together than 1%
            % as possible repeat eigenvalues (indicating eigenplanes).
            repeat_evals_idx = zeros(1,num_comps);
            for i = 1:num_comps-1
                tolerance = evals(i) / 100;
                repeat_evals_idx(i+1) = evals(i+1) > evals(i)-tolerance;
            end
            repeat_evals_idx = logical(repeat_evals_idx);
            
            % Drop components belonging to repeating eigenvalues (second component
            % of the pair only)
            comps2drop = find(repeat_evals_idx);
            repeat_evals = evals(repeat_evals_idx); % replace logical index with actual eigenvalues
            for i = 1:sum(repeat_evals_idx)
                evals = [evals(1:comps2drop(i)-1); evals(comps2drop(i)+1:end)];
                evecs = [evecs(:,1:comps2drop(i)-1), evecs(:,comps2drop(i)+1:end)];
                evecs_rel = [evecs_rel(:,1:comps2drop(i)-1), evecs_rel(:,comps2drop(i)+1:end)];
            end
            num_comps = num_comps - sum(repeat_evals_idx);
            
            %% Construct time series for components above significance threshold
            compts = zeros(num_comps, EEG.pnts, EEG.trials);
            
            % Note: component time series are scaled by relative eigenvalue,
            % facilitating amplitude comparisons within and between subjects.
            for c = 1:num_comps
                compts_temp = evecs_rel(:,c)' * ...
                    reshape(EEG.data, EEG.nbchan, []); % apply sensor weightings to concatenated trials
                compts(c,:,:) = reshape(compts_temp, ...
                    EEG.pnts, EEG.trials); % reshape back to individual trials
            end
            
            
            %% Save GED results to file
            
            % First, trim data from EEG, MEG, MEEG structs, to prevent saving and
            % loading from taking impractically long. The other analysis files tend
            % to primarily use the compts and the other struct fields anyway. Full
            % data sets can be loaded in when needed by calling loadMEEG with the
            % 'data' flag.
            EEG = rmfield(EEG, 'data');
            
            % Second, put results from GED in struct, to distinguish between
            % variables intrinsic to data and variables following from GED.
            GED.orig_evals = orig_evals; % for plotting
            GED.orig_evecs = orig_evecs; % for plotting
            GED.evals = evals;
            GED.repeat_evals = repeat_evals;
            GED.repeat_evals_idx = repeat_evals_idx;
            GED.evecs = evecs;
            GED.num_comps = num_comps;
            GED.compts = compts;
            GED.num_perm = num_perm;
            GED.perm_evals = obj.permresults.perm_evals;
            GED.covS = covS; % for plotting
            GED.covR = covR; % for plotting
            GED.kernel = obj.lastFilterKernal;
            GED.lo = obj.lastFilterLowerBound;
            GED.hi = obj.lastFilterUpperBound;
            
            disp('Saving results to file...');
            outputpath = obj.getPermutationDirectory();
            eegid = EEG.subject;
            GED_filename = fullfile(outputpath,[eegid '_GED.mat']);
            allrtidx = []; trialtype=[];
            save(GED_filename, 'EEG', 'allrtidx', ...
                'trialtype', 'GED');
            
            obj.lastGedResults = GED;
            
        end
        function obj = identifyMidFrontalComps( obj )
            
            if nargin < 2
                gedResults = load(obj.getGedResFileName);
                EEG = gedResults.EEG;
                GED = gedResults.GED;
            end
            
            midf = struct();
            
            %% Construct midfrontal theta template: Gaussian centered on FCz
            % Needs to be inside subject loop, because subject EEG.chanlocs
            % sometimes vary.
            % hard coded electrode
            
            midf.r2_cutoff = 0.5;
            
            fczidx = strcmpi('E6',{EEG.chanlocs.labels}); % E6 = FCZ
            eucdist = zeros(1,EEG.nbchan);
            
            for chani = 1:EEG.nbchan
                eucdist(chani) = sqrt( ...
                    (EEG.chanlocs(chani).X-EEG.chanlocs(fczidx).X)^2 + ...
                    (EEG.chanlocs(chani).Y-EEG.chanlocs(fczidx).Y)^2 + ...
                    (EEG.chanlocs(chani).Z-EEG.chanlocs(fczidx).Z)^2 );
            end
            
            midf.template = exp(-(eucdist.^2)/(2*50^2) );
            
            %% For each significant component, construct spatial forward filter model
            %  and calculate shared variance (r^2) with template
            
            % Create EEG topoplots + compute shared variance
            midf.template_r2 = zeros(1, GED.num_comps);
            midf.ffm_EEG = zeros(EEG.nbchan, GED.num_comps);
            for c = 1:GED.num_comps
                topo = GED.evecs(1:EEG.nbchan,c)' * ...
                    GED.covS(1:EEG.nbchan, 1:EEG.nbchan);
                
                % Flip component by sign of correlation with FCz template
                midf.ffm_EEG(:,c) = topo * sign(corr(topo', midf.template'));
                
                midf.template_r2(c) = corr(topo', midf.template')^2;
            end
            
            %% Select midfrontal components to use in following analyses
            comps = 1:GED.num_comps;
            midf.comps2use = comps(midf.template_r2 > midf.r2_cutoff);
            midf.num_comps = length(midf.comps2use);
            
            GED.midf = midf;
            obj.lastMidF = midf;
            obj.lastGedResults = GED;
            disp('Saving results to file...');
            outputpath = obj.getPermutationDirectory();
            eegid = EEG.subject;
            GED_filename = obj.getGedResFileName;
            allrtidx = []; trialtype=[];
            save(GED_filename, 'EEG', 'allrtidx', ...
                'trialtype', 'GED');
            
            
        end
        function mfcomps = getMidFrontalComps( obj )
            
            if nargin < 2
                gedResults = load(obj.getGedResFileName);
                EEG = gedResults.EEG;
                GED = gedResults.GED;
            end
            
            if isempty(obj.lastMidF)
                fprintf('Midfrontal Components Not Identified.');
                mfcomps = [];
            else
                mfcomps = GED.midf.comps2use;
            end
            
        end
        function obj = perform_gedperm_wrapper_parfor( obj, num_perm )
            % wrapper function for PERMUTEGED by  Marrit Zuure
            % available on github
            
            assert(obj.isSegmentedData, 'Requires segmented data.')
            % get filter kernal
            covS_window = obj.covS_window;
            covR_window = obj.covR_window;
            kernel = obj.lastFilterKernal;
            EEG = obj.getCurrentEEG;
            
            tic
            perm_evals = zeros(1,108,10);
            perm_settings = zeros(1,3005,10);
            perm_evecs = zeros(10,108,108);
            perm_covS = zeros(10,108,108);
            
            parfor pari = 1 : 1000 % 308.7 for 1000 permutations
                [perm_evals(:,:,pari),...
                    perm_settings(:,:,pari),...
                    perm_evecs(pari,:,:),...
                    perm_covS(pari,:,:)] = ...
                    permuteGED(EEG, covS_window, covR_window, kernel, 1, '95');
            end
            
            if strcmpi(method, '95')
                temp_evals = zeros(1,length(covSperm));
                for i = 1:length(covSperm)
                    temp_evals(i) = prctile(perm_evals(i,:), 95); % for each eigenvalue, extract 95th percentile from full set of permutations
                end
                perm_evals = temp_evals;
            end
            
            toc
            tic
            [perm_evals, perm_settings, perm_evecs, perm_covS] = ...
                permuteGED(EEG, covS_window, covR_window, kernel, 100, '95');
            toc
        end
        function [obj, perm_evals, perm_settings, perm_evecs, perm_covS] = ...
                perform_gedperm_wrapper( obj, num_perm )
            % wrapper function for PERMUTEGED by  Marrit Zuure
            % available on github
            assert(obj.isSegmentedData, 'Requires segmented data.')
            % get filter kernal
            covS_window = obj.covS_window;
            covR_window = obj.covR_window;
            kernel = obj.lastFilterKernal;
            EEG = obj.getCurrentEEG;
            
            eegid = EEG.subject;
            
            outputpath = obj.getPermutationDirectory();
            savefile = fullfile(outputpath,[eegid '_permtest.mat']);
            
            %% Check for existing permutation testing data; if not, run permutation test
            permtest_filename = savefile;
            if exist(permtest_filename,'file')
                load(permtest_filename, 'perm_settings');
                % Check whether loaded permutation data were generated using the
                % same settings (filter, window) as the GED. If not, the Z-scoring
                % and thus the significance thresholding will not be valid.
                % Obviously, MEEG also needs to be the same, but we don't check for
                % that.
                assert(all(perm_settings == [covS_window, covR_window, kernel]), ...
                    [ permtest_filename ' may not have been generated using the same ' ...
                    'filter kernel and time windows as were used for the GED.']);
                load(permtest_filename);
            else
                [perm_evals, perm_settings, perm_evecs, perm_covS] = permuteGED(EEG, covS_window, covR_window, kernel, num_perm, '95');
                save(permtest_filename, 'perm_evals', 'perm_settings', 'perm_evecs', 'perm_covS');
            end
            
            obj.permresults.perm_evals = perm_evals;
            
        end
        
        function savefile = getPermResFileName( obj )
            EEG = obj.getCurrentEEG;
            eegid = EEG.subject;
            outputpath = obj.getPermutationDirectory();
            savefile = fullfile(outputpath,[eegid '_permtest.mat']);
        end
        function savefile = getGedResFileName( obj )
            EEG = obj.getCurrentEEG;
            if ~isempty( EEG.subject )
                eegid = EEG.subject;
            else
                eegid = obj.s.subj_basename;
            end
            outputpath = obj.getPermutationDirectory();
            savefile = fullfile(outputpath,[eegid '_GED.mat']);
        end
        function compts = getCompTimeSeries( obj, compidx )
            gedResults = load(obj.getGedResFileName);
            GED = gedResults.GED;
            EEG = obj.getCurrentEEG;
            if nargin < 2 || isempty(compidx)
                compidx = [1:size(GED.compts,1)];
            end
            
            compts = GED.compts(compidx,:,:);
            
        end
        function EEG = getCompTimeSeriesSetFile( obj, compidx )
            if nargin < 2
                compidx = [];
            else
            end
            EEG = obj.getCurrentEEG;
            compts = obj.getCompTimeSeries(compidx);
            EEG.data = compts;
            EEG = eeg_checkset(EEG);
            EEG = eeg_checkchanlocs(EEG);
            for chani=1:length(compidx)
                EEG.chanlocs(chani).labels = compidx(chani);
                EEG.chanlocs(chani).type = 'EEG';
                EEG.chanlocs(chani).ref = 'GEDcomp';
            end
            
        end
        function obj = ged( obj )
            % perform checks (exists, same size)
            % [evecs, evals] = eig(double(obj.SCovMat), double(obj.RCovMat));
            MEEG = obj.s.EEG;
            %% Construct frequency filter kernel
            nyquist = MEEG.srate/2;
            
            % GED filter settings
            filter_lo = 5; % lower theta bound in Hz
            filter_hi = 7; % upper theta bound in Hz
            filter_order = 3000; % higher order is better frequency resolution, poorer temporal resolution
            
            % GED time window (in ms) to use in constructing covariance matrices
            covS_window = [-1000 1000]; % signal matrix
            covR_window = [-1000 1000]; % noise matrix
            kernel = fir1(filter_order, [filter_lo/nyquist filter_hi/nyquist]);
            obj.lastFilterKernal = kernal;
            
            
            %% Perform a generalized eigendecomposition on theta-filtered vs. unfiltered data
            [evals, evecs, evecs_rel, covS, covR, evecs_unnormed] = ...
                performGED(MEEG, covS_window, covR_window, kernel);
            num_perm = 50;
            
            %             MEEG.data = gpuArray(MEEG.data);
            %             MEEG.data = gather(MEEG.data);
            tic;
            [perm_evals, perm_settings, perm_evecs, perm_covS] = ...
                permuteGED(MEEG, covS_window, covR_window, kernel, num_perm, '95');
            toc;
            save(permtest_filename, 'perm_evals', 'perm_settings', 'perm_evecs', 'perm_covS');
            toc;
            %% Determine significant components
            num_comps = sum(evals > perm_evals');
            %% Cleanup: Drop components that potentially represent eigenplanes
            
            % Flag significant eigenvalues that are closer together than 1%
            % as possible repeat eigenvalues (indicating eigenplanes).
            repeat_evals_idx = zeros(1,num_comps);
            for i = 1:num_comps-1
                tolerance = evals(i) / 100;
                repeat_evals_idx(i+1) = evals(i+1) > evals(i)-tolerance;
            end
            repeat_evals_idx = logical(repeat_evals_idx);
            
            % Drop components belonging to repeating eigenvalues (second component
            % of the pair only)
            comps2drop = find(repeat_evals_idx);
            repeat_evals = evals(repeat_evals_idx); % replace logical index with actual eigenvalues
            for i = 1:sum(repeat_evals_idx)
                evals = [evals(1:comps2drop(i)-1); evals(comps2drop(i)+1:end)];
                evecs = [evecs(:,1:comps2drop(i)-1), evecs(:,comps2drop(i)+1:end)];
                evecs_rel = [evecs_rel(:,1:comps2drop(i)-1), evecs_rel(:,comps2drop(i)+1:end)];
            end
            num_comps = num_comps - sum(repeat_evals_idx);
            
            %% Construct time series for components above significance threshold
            compts = zeros(num_comps, MEEG.pnts, MEEG.trials);
            
            % Note: component time series are scaled by relative eigenvalue,
            % facilitating amplitude comparisons within and between subjects.
            for c = 1:num_comps
                compts_temp = evecs_rel(:,c)' * reshape(MEEG.data, MEEG.nbchan, []); % apply sensor weightings to concatenated trials
                compts(c,:,:) = reshape(compts_temp, MEEG.pnts, MEEG.trials); % reshape back to individual trials
            end
            
            
            figure; plot(reshape(squeeze(compts(1,:,:)), 1, []))
            
            orig_evals = evals;
            orig_evecs = evecs_unnormed;
            
            
            % sort so largest component is first
            [evals,sidx]  = sort( diag(evals),'descend' );
            evecs = evecs(:,sidx);
            
            obj.lastEvecs = evecs;
            obj.lastEvals = evals;
            figure;
            subplot(211)
            obj.plotComponentLambda( evals, evecs_rel);
            subplot(212)
            obj.plotComponentLambda( perm_evals, perm_evecs);
            
            
        end
        function obj = getEigComponent( obj, compnumber )
            
            if nargin < 2
                overrideMaxComponent = false;
            else
                overrideMaxComponent = true;
                obj.SCovCompNumber = compnumber;
            end
            
            evals = obj.lastEvals;
            evecsT = obj.lastEvecs;
            thetacov = obj.SCovMat;
            thetafilt = obj.filtData;
            
            % find best component and compute filter projection
            [B,maxcomp] = sort(diag(evals));
            thetamap    = thetacov*evecsT;
            if overrideMaxComponent == true && compnumber < length(thetamap)
                thetamap    = thetamap(:,compnumber);
            else
                thetamap    = thetamap(:,maxcomp(1));
            end
            obj.maxcomp = maxcomp;
            % fix sign of map (max is positive)
            [~,maxe] = max(abs(thetamap));
            thetamap = thetamap * sign(thetamap(maxe));
            
            % theta time series component
            if overrideMaxComponent == true && compnumber < length(thetamap)
                compno_select = compnumber;
                evecsT_select = evecsT(:,compnumber);
            else
                compno_select = maxcomp(end);
                evecsT_select = evecsT(:,maxcomp(end));
            end
            if thetafilt.trials > 1
                thetacomp = reshape(thetafilt.data, ...
                    [thetafilt.nbchan thetafilt.trials *thetafilt.pnts])' ...
                    * evecsT_select;
                % bbdata = reshape(obj.s.EEG.data, ...
                %    [obj.s.EEG.nbchan obj.s.EEG.trials *obj.s.EEG.pnts])'
                times = 0:...
                    1/thetafilt.srate: ...
                    thetafilt.pnts * thetafilt.trials / ...
                    thetafilt.srate;
                obj.CompTimes = times;
            else
                thetacomp = thetafilt' * evecsT_select;
                obj.CompTimes = thetafilt.times;
                % obj.bbdata = bbdata;
            end
            
            % fix sign of time series according to sign of correlation with EEG
            thetacomp = thetacomp ...
                * sign(corr(thetacomp,...
                filterFGx(obj.s.EEG.data(maxe,:),...
                obj.s.EEG.srate,...
                obj.lastFilterHz,...
                obj.lastFilterFhwm)'));
            
            obj.SCovComp = thetacomp;
            
            
        end
        function obj = troughCovMat( obj )
            saveplot = obj.saveplot;
            
            % starting parameters
            EEG = epoch2continuous(obj.s.EEG);
            filtfreq   = obj.lastFilterHz;
            fwhm        = obj.lastFilterFhwm;
            
            % get filtered cov component + details
            filtcomp.signal   = obj.SCovComp; % from create cov functions
            filtcomp.pnts     = numel(filtcomp.signal);
            filtcomp.times    = EEG.times;
            filtcomp.srate    = EEG.srate;
            filtcomp.ssignal_m  = smooth(filtcomp.signal);
            
            % pipeline to identify trough
            
            % determine window size around trough
            windowSize = 1/4; % window size is 1/4 cycle (1/8 of either side)
            nwin = ceil(filtcomp.srate / filtfreq * windowSize / 2);
            
            filtcomp.nwin = nwin;
            
            % find trough across entire filtered signal
            % algos compared: diff > findpeaks
            % performed on smoothed data
            
            % find large outliers that distrupt mean/sd calc
            % for power threshold
            
            filtcomp.ssignal_m_outlier = filloutliers(filtcomp.ssignal_m, ...
                'center','mean','ThresholdFactor', 3);
            
            % A covariance matrix of broadband data was constructed
            % around alpha peaks (only peaks corresponding to >1
            % standard deviation above the mean alpha peak amplitudes
            % were included; this limits the analysis to periods of
            % high alpha power)
            
            % calculate thresholds
            powerSdMultiplier = 1; % one SD
            filtcomp.meanpow = mean( filtcomp.ssignal_m_outlier .^2 );
            filtcomp.sdpow = std(  filtcomp.ssignal_m_outlier .^2 );
            filtcomp.threspow = filtcomp.meanpow + filtcomp.sdpow * powerSdMultiplier;
            filtcomp.thresamp = -1 * sqrt(abs(filtcomp.threspow)); % amplitude for plotting
            
            % find all troughs (inverse to find peaks)
            filtcomp.troughidx = ...
                find(diff(sign(diff( filtcomp.ssignal_m_outlier )))>0)+1;
            
            % isolate troughs to > 1 SD power
            filtcomp.troughidxThresPow = ...
                filtcomp.ssignal_m(filtcomp.troughidx).^2 ...
                > filtcomp.threspow;
            
            validTroughIdx =  ...
                filtcomp.troughidx(  filtcomp.troughidxThresPow );
            
            troughs = validTroughIdx;
            filtcomp.validTroughIdx = troughs;
            
            % remove start and end of data to prevent
            % errors on windowing
            
            troughs(troughs<nwin+1) = [];
            troughs(troughs>filtcomp.pnts-nwin-1) = [];
            
            % prepare matrix, must be channel covariance
            covT = zeros(EEG.nbchan);
            
            % trough-locked covariance
            for ti=1:length(troughs)
                tmpdat = EEG.data(:,troughs(ti)-nwin:troughs(ti)+nwin);
                tmpdat = bsxfun(@minus,tmpdat,mean(tmpdat,2));
                covT   = covT + (tmpdat*tmpdat')/nwin;
            end
            
            covA = cov(EEG.data');
            
            covT = covT./ti;
            
            covT = zeros(EEG.nbchan);
            
            % trough-locked covariance
            for ti=1:length(troughs)
                tmpdat = EEG.data(:,troughs(ti)-nwin:troughs(ti)+nwin);
                tmpdat = bsxfun(@minus,tmpdat,mean(tmpdat,2));
                covT   = covT + (tmpdat*tmpdat')/nwin;
            end
            covT = covT./ti;
            
            obj.pac_covT = covT; % trough broadband
            
            
            %% GED to get gamma peak/trough networks
            % examine difference betwee trough bb and all broadband
            % make sure it is double computation
            
            [evecs,evals] = eig(double(obj.pac_covT),double(obj.RCovMat));
            
            [evals,sidx]  = sort( diag(evals),'descend' );
            evecs = evecs(:,sidx);
            
            [~,compidx]   = sort(diag(evals)); % max component
            
            figure(1), clf
            subplot(211)
            plot(evals,'ks-','markersize',10,'markerfacecolor','w')
            xlabel('Component')
            ylabel('\lambda')
            
            plot(1:128,sort(diag(evals)))
            maps    = obj.pac_covT*evecs; % forward model of filter
            gamnet1 = maps(:,compidx(end));
            
            % fix sign
            hz = linspace(0,filtcomp.srate,numel(gamnet1));
            
            [~,idx] = max(abs(gamnet1));
            gamnet1 = gamnet1 * sign(gamnet1(idx));
            % plot channel power spectra (hard-coded for POz)
            subplot(332)
            figure;
            
            
            
            for xi = 1 : 20
                subplot(5,4,xi)
                plot(hz,abs(fft(maps(:,compidx(end-xi-1))))/numel(gamnet1),'k')
                set(gca,'xlim',[0 90],'xtick',10:20:90), axis square
                xlabel('Frequency (Hz)'), ylabel('Power (a.u.)')
                title('Trough Component')
            end
            plot(filtcomp.times,maps(:,compidx(end)),'k')
            figure;
            for xi = 1 : 20
                subplot(5,4,xi)
                plot(filtcomp.times,compidx(end-xi-1),'k')
                set(gca,'xlim',[10 100]), axis square
                xlabel('Frequency (Hz)'), ylabel('Power (a.u.)')
                title('Trough Component')
            end
            
            plot(gamnet1)
            
            obj.SCovMat = gamnet1;
            
            if saveplot == true
                % debug code to visualize peak detection performance
                figure(34);
                subplot(2,2,1:2);
                % original signal
                plot(filtcomp.times, filtcomp.signal, 'y', 'LineWidth', 3 )
                xlim([10*filtcomp.srate 50*filtcomp.srate])
                hold on;
                
                % plot threshold line (SD)
                yline(filtcomp.thresamp, 'r');
                % smoothed signal
                plot(filtcomp.times, filtcomp.ssignal_m_outlier);
                % all troughs
                plot(filtcomp.times(filtcomp.troughidx), ...
                    filtcomp.signal(filtcomp.troughidx),'g*');
                % troughs > threshold
                times_t = filtcomp.times( validTroughIdx );
                amp_t = filtcomp.signal( validTroughIdx );
                plot(times_t, ...
                    amp_t,...
                    'ko','MarkerSize',10);
                subplot(2,2,3);
                imagesc(covA), axis square;
                title('Broadband Cov');
                
                subplot(2,2,4);
                imagesc(covT), axis square;
                title('Broadband Cov tied to Troughs');
            end
            
            
        end
        function filt_order = createFiltOrder(obj, lowerBound, upperBound)
            % inputs
            TRANSWIDTHRATIO = 0.25;
            filter_lo        = lowerBound; % lower band Hz
            filter_hi        = upperBound; % higher bands Hz
            edgeArray = sort([filter_lo filter_hi]);
            nq = obj.currentSampleRate / 2;
            
            maxTBWArray = edgeArray; % Band-/highpass
            maxTBWArray(end) = nq - edgeArray(end);
            maxDf = min(maxTBWArray);
            df = min([max([edgeArray(1) * TRANSWIDTHRATIO 2]) maxDf]);
            
            if ~exist('filtOrder','var')
                df = min([max([edgeArray(1) * TRANSWIDTHRATIO 2]) maxDf]);
                filtOrder = 3.3 / (df / obj.currentSampleRate); % Hamming window
                filtOrder = ceil(filtOrder / 2) * 2; % Filter order must be even.
                filt_order = filtOrder*3;
            else
            end
        end
        function kernal = constructFilterKernal(obj, lowerBound, upperBound, filterOrder)
            if nargin < 2
                % GED filter settings
                lowerBound = 5; % lower theta bound in Hz
                upperBound = 7; % upper theta bound in Hz
                filter_order = 3000; % higher order is better frequency resolution, poorer temporal resolution
            end
            
            nq = obj.currentSampleRate / 2;
            filter_order = obj.createFiltOrder(lowerBound, upperBound);
            
            % calculate filt order for bandpass only
            % adapted from pop_eegfiltnew (EEGLAB)
            
            %filter_order = 3000;
            kernal = fir1(filter_order, [lowerBound/nq upperBound/nq]);
            obj.lastFilterKernal = kernal;
            obj.lastFilterOrder = filter_order;
            obj.lastFilterLowerBound = lowerBound;
            obj.lastFilterUpperBound = upperBound;
            
        end
        function obj = setGedTimeWindows(obj, covS_window, covR_window)
            assert(isvector(covS_window) && isvector(covR_window));
            % GED time window (in ms) to use in constructing covariance matrices
            obj.covS_window = covS_window*1000; % signal matrix
            obj.covR_window = covR_window*1000; % noise matrix
            
            fprintf('cov S window set: %s \ncov R window set: %s', ...
                num2str(covS_window), num2str(covR_window));
        end
    end
    methods % statistical methods
    end
    methods % GED Bounds
        function vec = gausFxSmooth(obj, peakf, hz, fwhm )
            peakf=2
            fwhm=5;
            if peakf == hz(1)
                peakf = hz(450);
            end
            impulsevec = zeros(1, length(hz));
            impulsevec(hz==peakf)=1;
            
            
            % frequency-domain Gaussian
            s  = fwhm*(2*pi-1)/(4*pi); % normalized width
            x  = hz-peakf;             % shifted frequencies
            fx = exp(-.5*(x/s).^2);    % gaussian
            
            dataX = impulsevec;
            % IFFT
            figure; plot(hz, abs( dataX.*fx ));
            xlim([50 60]);
            %     hold on;
            %     xline(peakf);
            %     xlim([25 30])
            %     convres = 2*real( ifft( dataX.*fx ));
            %
            %     figure; plot(fx);
            %     figure; plot(hz,convres);
            
            
        end
        function obj = gedBounds( obj, freqRange, numfrex )
            % computation outline
            % define parameters
            % retreive data
            % transform data (example)
            % results: preallocation
            % results: computations
            % results: transform results (examples)
            % store results
            
            % function defaults
            if nargin < 2, lowfreq = 2; highfreq = 80; numfrex = 200;
            else, lowfreq = freqRange(1); highfreq = freqRange(2); 
            end
            
            % 01 Dataset Import
            % dataset is loaded in EEG property of this object
            
            % 02 Dataset Transform
            
            EEG = obj.getCurrentEEG;    % retrieve EEG set
            preEEG.trials = EEG.trials;
            preEEG.pnts = EEG.pnts;
            obj.epoch2continuous;       % conv. to continuous
            EEG = obj.getCurrentEEG;    % retrieve EEG set
            EEG.data = double(EEG.data);
            % onset times for epoching resting-state data
            % start at 2s, step by 2s, end at last poin minus 4 s
            onsets = EEG.srate*2:2*EEG.srate:EEG.pnts-EEG.srate*4;
            snipn  = 2*EEG.srate;  % note edges
            
            % 03a Results: Parameters
            filterType = 'fgx'; % fgx or filtfilt or eegfilt

            stds = linspace(2,5,numfrex);       % std for gaussian
            
            frex = logspace(log10(lowfreq),...  % freq. vector
                log10(highfreq),numfrex);
%             figure; plot(stds);
%             figure; plot(frex);
            
            % 03b Results: Preallocation            
            res = struct();
            [evals,evecs,maps] = deal(zeros(numfrex,EEG.nbchan));
            R = zeros(length(onsets),EEG.nbchan,EEG.nbchan);
            
            % 03c Results: Computation
            for segi=1:length(onsets) % create & clean  R covariance matrix
                snipdat = EEG.data(:,onsets(segi):onsets(segi)+snipn);
                snipdat = bsxfun(@minus,snipdat,mean(snipdat,2));
                R(segi,:,:) = snipdat*snipdat'/snipn;
            end
            % clean R
            meanR = squeeze(mean(R));
            dists = zeros(1,size(R,1));
            for segi=1:size(R,1)
                r = R(segi,:,:);
                dists(segi) = sqrt( sum((r(:)-meanR(:)).^2) );
            end
            R = squeeze(mean( R(zscore(dists)<3,:,:) ,1));
            % regularized R
            gamma = .01;
            Rr = R*(1-gamma) + eye(EEG.nbchan)*gamma*mean(eig(double(R)));

            %% create S narrowband matrices
            
            for fi=1:numfrex
                
               % if frex(fi) > 10
                   % fprintf('target freq');
              %  end
                % filter data
                %obj.filtfilt(frex(fi),stds(fi)/2, EEG);
                %fdat = obj.filtData.data;
               %tic;
               fprintf('s: %s f: %2.1f s: %2.1f\n',  obj.s.subj_basename, frex(fi),std(fi));
                 [fdat empVals] = ...
                     filterFGx(EEG.data,EEG.srate,...
                     frex(fi),stds(fi));
               %toc;
                 % figure; plot(fdat(1, 1:1000))
               % confirm filtered data peak at target
               %                    powrw = pwelch(fdat(34,:), 500);
               %                 figure; plot(mean(powrw,2)); xlim([0,15]);
               
%                 fdat2 = reshape(permute(fdat, [2 1]), [], preEEG.trials , EEG.nbchan);
%                 fdat2 = permute(fdat2, [3,2,1]);
%                % figure; plot(squeeze(fdat2(34,1,1:1000)));                                                
%                [pxx,f] = pwelch(squeeze(fdat2(34,1,:)),500,300,500,EEG.srate);
%                [pxx,f] = pwelch(squeeze(fdat(34,:)),500,300,500,EEG.srate);
% [pxx,w] = pwelch(squeeze(fdat2(34,1,:)),[],[],500);
% [powr2, snr2, hz] = gedTsPow( squeeze(fdat2(34,:,:)), EEG.srate, ...
%                     preEEG.pnts, [0 20], frex(fi) );
% [powr2, snr2, hz] = gedTsPow( squeeze(fdat(34,:)), EEG.srate, ...
%                     preEEG.pnts, [0 20], frex(fi) );
%                 %figure; plot(fdat);
                
                % calculate power spectrum
                
                
                % preallocation of S
                S = zeros(length(onsets),EEG.nbchan,EEG.nbchan);
                
                % full S cov                
                for segi=1:length(onsets)
                    snipdat = fdat(:,onsets(segi):onsets(segi)+snipn);
                    snipdat = bsxfun(@minus,snipdat,mean(snipdat,2));
                    S(segi,:,:) = snipdat*snipdat'/snipn;
                end
                
                % clean S
                meanS = squeeze(mean(S));
                dists = zeros(1,size(S,1));
                for segi=1:size(S,1)
                    se = S(segi,:,:);
                    dists(segi) = sqrt( sum((se(:)-meanS(:)).^2) );
                end
                S = squeeze(mean( S(zscore(dists)<3,:,:) ,1));
                
                % global variance normalize
                S = S / (std(S(:))/std(R(:)));

                % check covar mat via plots
              %  obj.plotEigMatrices(S, R, frex(fi), stds(fi)/2);
                % compute GED
                % normalize covariance matrix between 0-1
                
                [W,L] = eig(S,Rr);
                [evals(fi,:),sidx] = sort(diag(L),'descend');
                W = W(:,sidx);
                
                
                
                %% Normalize eigenvectors to unit length
                % Also create variable containing eigenvectors normalized to the relative
                % eigenvalue. This is later used to scale the component time series, to
                % facilitate power and amplitude comparisons between components and
                % between subjects.
                evecs_normed = W;
                evecs_unnormed = evecs;
                evecs_rel = zeros(size(evecs_normed));
                for v = 1:size(evecs_normed,2)
                    evecs_normed(:,v) = evecs_normed(:,v)/norm(evecs_normed(:,v)); % normalize to unit length
                    rel_eval = evals(fi,v)/sum(evals(fi)); % extract relative eigenvalue
                    evecs_rel(:,v) = evecs_normed(:,v) * rel_eval; % normalize to relative eigenvalue
                    evals_rel(:,v) = rel_eval;
                end
                %imagesc(evecs_rel); axis square;
                res.evecsdb.evals_rel{fi} = evals_rel;
                res.evecsdb.evecs_rel{fi} = evecs_rel;
                res.evecsdb.evecs_unnormed{fi} = evecs_unnormed;
                
                % store top component map and eigenvector
                maps(fi,:) = W(:,1)'*S;
                evecs(fi,:) = W(:,1);
                
                %% Construct time series for components above significance threshold
                % compts = zeros(num_comps, EEG.pnts);
                res.compts(fi, :) = ...
                    evecs_rel(:,1)' * EEG.data;

              %   powrw = pwelch(res.compts(fi, :), 500); 
                %figure; plot(mean(powrw,2)); xlim([0,15]);
%                 figure; plot(fdat(); xlim([0,3000]);
            end
            
            %% correlation matrices for clustering
            
            E = zscore(evecs,[],2);
            evecCorMat = (E*E'/(EEG.nbchan-1)).^2;
            
            res.maps = maps;
            res.evecs = evecs;
            res.evals = evals;
            res.evecsCorMatZ = evecCorMat;
            
            
            % apply sensor weightings to concatenated trials
            % Note: component time series are scaled by relative eigenvalue,
            % facilitating amplitude comparisons within and between subjects.
            
            
            
            %% determine the optimal epsilon value
            
            % range of epsilon parameter values
            nepsis = 50;
            epsis  = linspace(.001,.05,nepsis);
            qvec   = nan(nepsis,1);
            
            for epi=1:length(epsis)
                
                % scan
                freqbands = dbscan(evecCorMat, epsis(epi),3,'Distance','Correlation');
                if max(freqbands)<4, continue; end
                
                % compute q
                qtmp = zeros(max(freqbands),1);
                MA = false(size(evecCorMat));
                for i=1:max(freqbands)
                    M = false(size(evecCorMat));
                    M(freqbands==i,freqbands==i) = 1;
                    qtmp(i) = mean(mean(evecCorMat(M))) / ...
                        mean(mean(evecCorMat(~M)));
                    MA = MA+M;
                end
                qvec(epi) = mean(qtmp) + log(mean(MA(:)));
            end
            
            res.epsis = epsis;
            res.qvec = qvec;
            
            % run it again on the best epsilon value
            [~,epsiidx] = findpeaks(qvec,'NPeaks',1,'SortStr','descend');
            if isempty(epsiidx), epsiidx = round(nepsis/2); end
            freqbands = dbscan(evecCorMat,epsis(epsiidx),3,'Distance','Correlation');
            
            res.bestepsis = epsis(epsiidx);
            res.bestepsisIdx = epsiidx;
            res.bestqvec = qvec(epsiidx);
            
            % dissolve tiny clusters, and renumber all clusters consecutively
            newc = cell(4,1); n=1;
            for i=1:max(freqbands)
                cc = bwconncomp(freqbands==i);
                for ci=1:cc.NumObjects
                    if length(cc.PixelIdxList{ci})>2
                        newc{n} = cc.PixelIdxList{ci};
                        n = n+1;
                    end
                end
            end
            
            % occasionally large clusters will be out of order
            % sorted by first index of each cluster  
            for fix_i =  1 : length(newc)
                tmpcell = newc{fix_i};
                startvec(fix_i,1) =fix_i;  % original index
                startvec(fix_i,2) =tmpcell(1); % first element of each cluster
            end
            newsort = sortrows(startvec, 2);
            newc = newc(newsort(:,1));  % use original index to reorder newc
            
            
            freqbands = -ones(size(frex));
            for ni=1:n-1
                freqbands(newc{ni}) = ni;
            end
            
            % identify lower and upper bounds of identified frequency bands
            
            res.freqbands = freqbands;
            res.nfreqbands = max(freqbands);
            res.frex = frex;
            
            chanlocs = EEG.chanlocs;
            res.chanlocs = chanlocs;
            res.evecCorMat = evecCorMat;
            %% average correlation coefficient within each cluster
            
            avecorcoef = zeros(max(freqbands),2);
            for i=1:max(freqbands)
                submat = evecCorMat(freqbands==i,freqbands==i);
                avecorcoef(i,1) = mean(nonzeros(tril(submat,-1)));
                avecorcoef(i,2) = mean(frex(freqbands==i));
            end
            
            res.avecorcoef = avecorcoef;
            res.filterType = filterType;
            obj.gedBoundsRes = res;
            
            showplots = false;
            if showplots
                obj.plotBandBoundries;
                obj.plotPcaTopography;
            end
            
            
            
            % 04 Export Results
            res.srate = EEG.srate;
            res.duration = EEG.xmax - EEG.xmin;
            
        end
    end    
    methods % visualization methods
        function obj = plotEigMatrices(obj, S, R, frex, stds)
            figure;
            chanAxisLabel = 'Channel Number';
            subplot(1,2,1)
            imagesc(zscore(S), [-4 4])
            title(sprintf('S Matrix (%d +/- %dHz)', frex, stds))
            axis square
            xlabel(chanAxisLabel)
            ylabel(chanAxisLabel)
            subplot(1,2,2)
            imagesc(zscore(R), [-4 4])
            title('R Matrix')
            axis square
            xlabel(chanAxisLabel)
            ylabel(chanAxisLabel)
            colormap(jet)
            sgtitle(sprintf('Narrow-Band & Broadband Covariance Matrices\n%s', obj.s.subj_basename), 'Interpreter', 'none');
            set(gcf, 'Position', [680,687,560,291])
            
               fname = fullfile(obj.pathdb.figures, ...
                   sprintf('EigMatrix_%s_%d_HZ.png', obj.s.subj_basename, frex(fi)));                
            obj.filedb.eigplot = fname;
            saveas(gcf,fname);
           % close gcf
            
        end
        function obj = plotDataTimeSeries(obj)
            EEGBrowser(obj.s.EEG);
        end
        function obj = plotLastFiltTimeSeries(obj)
            EEGBrowser(obj.filtData);
        end
        function fig = plotCovMat( obj, covmat )
            % use last cov mat if not available
            if nargin < 2
                covmat = obj.lastCovMat;
            end
            figure(1);
            imagesc(covmat), axis square;
            fig = figure(1);
        end
        function fig = plotLastCovMatOutliers( obj )
            EEG = obj.s.EEG;
            data2 = obj.lastEEG.data;
            
            covdistz = obj.lastCovdistz;
            % show the covariance distances
            figure(1), clf
            subplot(2,3,1:2)
            plot(covdistz,'ks-','linew',2,'markerfacecolor','w','markersize',12)
            xlabel('Trial'), ylabel('Z_{dist}')
            title('Z-scored covariance distances')
            
            
            % histogram of distances
            subplot(233)
            histogram(covdistz,10)
            xlabel('Distances'), ylabel('Count')
            title('Histogram of distances');
            
            % plot time courses
            subplot(212), hold on
            plot(EEG.times,mean(EEG.data(31,:,:),3),'k','linew',2)
            plot(EEG.times,mean(data2(31,:,:),3),'r','linew',2)
            
            % make the plot look a bit nicer
            xlabel('Time (a.u.)')
            legend({'Original data';'Trials removed'})
            title('Time series before and after covariance cleaning')
            zoom on
            
            fig = figure(1);
        end
        function fig = plotLastEigComponent( obj, compnumber )
            if nargin < 2
                compnumber = obj.maxcomp(end);
            else
                obj.SCovCompNumber = compnumber;
                fprintf("\n\nSwitch to component #%d.\n\n", compnumber);
            end
            
            obj.getEigComponent( compnumber );
            
            zoomtime = obj.s.EEG.srate * 10; % 10 seconds
            zoomstart = randi([zoomtime ...
                size(obj.CompTimes,2)]);
            zoomindex = zoomstart:zoomstart+zoomtime;
            
            figure(1);
            subplot(2,1,1)
            plot(obj.CompTimes(1:end-1), obj.SCovComp);
            title(sprintf('Filtered Component %d',...
                compnumber))
            subplot(2,1,2)
            plot(obj.CompTimes(zoomindex), ...
                obj.SCovComp(zoomindex,:));
            title(sprintf('10 second zoom view'));
            fig = figure(1);
            
        end
        function fig = plotComponentLambda( obj, evals, evecs)
            % initialize maps and time series
            if nargin < 2
                evecs = obj.lastEvecs;
                evals = obj.lastEvals;
            end
            EEG = obj.s.EEG;
            
            maps   = zeros(4,EEG.nbchan);
            compts = zeros(4,EEG.pnts,EEG.trials);
            figure(1), clf
            subplot(211);
            plot(evals,'ks-','markersize',10,'markerfacecolor','w');
            xlabel('Component Number');
            ylabel('\lambda');
            
            
            for compi=1:4
                
                % compute map and adjust sign
                maps(compi,:) = evecs(:,compi)'*obj.SCovMat;
                [~,idx]       = max( abs(maps(compi,:)) );
                maps(compi,:) = maps(compi,:) * sign(maps(compi,idx));
                
                % compute time series
                tmpts = evecs(:,compi)'*reshape(EEG.data,EEG.nbchan,[]);
                compts(compi,:,:) = reshape( tmpts,[EEG.pnts EEG.trials] );
                
                subplot(2,4,4+compi)
                topoplotIndie(maps(compi,:),EEG.chanlocs,'electrodes','off','numcontour',0);
                title([ 'Component ' num2str(compi) ])
            end
            fig = figure(1);
            % debug
            %             [evals,sidx]  = sort( diag(evals), 'descend' );
            %             evecs = evecs(:,sidx);
            
            
            
            
            
            
            %[~,compidx]   = sort(diag(evals)); % max component
        end
        function fig = createScreePlot( obj, GED )
            
            if nargin < 2
                
                gedResults = load(obj.getGedResFileName);
                GED=gedResults.GED;
                EEG = gedResults.EEG;
            end
            
            % create Scree Plot
            evals_to_plot = GED.evals(:);
            
            comps = 1:length(evals_to_plot);
            perm_evals = evals_to_plot(evals_to_plot > ...
                GED.perm_evals');
            mfcomps = GED.evals(obj.getMidFrontalComps)';
            titlestr = sprintf('%s (%s): Sig. Comps %d/%d Perm: %d', ...
                gedResults.EEG.subject, ...
                gedResults.EEG.group, ...
                length(evals_to_plot),...
                GED.num_comps, ...
                GED.num_perm);
            
            figure(1);
            plot(comps, evals_to_plot, 'b.', 'MarkerSize', 25);
            hold on;
            plot( comps(evals_to_plot > ...
                GED.perm_evals'), perm_evals, ...
                'ro',  'MarkerSize', 15);
            plot( comps(obj.getMidFrontalComps), mfcomps , ...
                'g.',  'MarkerSize', 15);
            
            plot(comps,  GED.perm_evals', '--k');
            % legend('eigenvalues', 'significance threshold');
            xlabel('Components');
            ylabel('S to R power ratio (\lambda)');
            ylim([-.15 max(evals_to_plot)+0.1]);
            xlim([-5 length(evals_to_plot)]); % hard-coding number of sensors
            set(gca, 'Box', 'on');
            
            ax = gca;
            ax.FontSize = 20;
            xticks('auto');
            yticks('auto');
            title(titlestr, 'Interpreter', 'none');
        end
        function fig = createCompSpectrogram( obj, compidx )
            gedResults = load(obj.getGedResFileName);
            GED = gedResults.GED;
            if nargin < 2 || isempty(compidx)
                compidx = [1:size(GED.compts,1)];
            end
            pop_spectopo(obj.getCompTimeSeriesSetFile(compidx));
        end
        function res = getGedBoundsRawVec( obj, showplots )
            if nargin < 2
                showplots = false;
            end
            
            gRes = obj.gedBoundsRes;
            frex = gRes.frex;
            freqbands = gRes.freqbands;
            
            % result variables
            res = {};
            tbnds_all = {};
            template_vec = zeros(1, size(frex,2));
            ivector=[];
            [ivector.range, ...
                ivector.edge, ...
                ivector.lower, ...
                ivector.upper, ...
                ivector.loclust, ...
                ivector.upclust, ...
                ivector.lofrex, ...
                ivector.upfrex] = deal(template_vec);
            
            for i=1:max(freqbands)
                tbnds = frex(freqbands==i);
                tbnds = dsearchn(frex',tbnds([1 end])');
                ivector.range(tbnds(1):tbnds(2)) = 1;
                ivector.edge([tbnds(1) tbnds(2)]) = 1;
                ivector.lower(tbnds(1)) = 1;
                ivector.upper(tbnds(2)) = 1;
                ivector.loclust(tbnds(1)) = i; % assign cluster
                ivector.upclust(tbnds(2)) = i;
                ivector.lofrex(tbnds(1)) = frex(tbnds(1));
                ivector.upfrex(tbnds(2)) = frex(tbnds(2));
            end
            
            obj.gedBoundsRes.ivector = ivector;
            res = ivector;
            
            publication_plot = false;
            if publication_plot
                f= figure; 
                title('Range of Empirical Frequency Boundries');
                xlabel('Frequency (hz)');
                l = frex(find(ivector.lower));
                u = frex(find(ivector.upper));                
                c = gRes.avecorcoef(:,1)';
                
                arrayfun(@(a)xline(a, 'LineWidth', 1, ...
                    'Color', 'c'), l); 
                arrayfun(@(a)xline(a, 'LineWidth', ...
                    1, 'Color', 'b'), u); 
                arrayfun(@(a)...
                    patch([l(a) l(a) u(a) u(a)],...
                    [0 c(a) c(a) 0], c(a),'EdgeColor','none', ...
                        'LineWidth',1,'FaceAlpha',1), ...
                        [1:length(l)])         
                   colorbar; 
                   colormap(flipud(bone));                  
                   
                   caxis([.8 1]);
                   xlim([0 80]);
                   ylim([.5 1]);
                   pbaspect([4 1 1])
                   
            end
            
            if showplots
                %figure(1);
                figure;
                subplot(4,1,1); plot(frex, ivector.range);
                title('Range of Bounds');
                subplot(4,1,2); plot(frex, ivector.edge);
                title('Edges of Bounds');
                subplot(4,1,3); 
                %plot(frex, ivector.lower,'r'); hold on;
                %plot(frex, ivector.upper,'b');
                arrayfun(@(a)xline(a, 'LineWidth', 1, 'Color', 'c'), frex(find(ivector.lower))); 
                arrayfun(@(a)xline(a, 'LineWidth', 1, 'Color', 'b'), frex(find(ivector.upper))); 

                l = frex(find(ivector.lower))
                u = frex(find(ivector.upper))
                [l; u]'
%                 x = [l(9) l(9) u(9) u(9) ];
%                 y = [0 1 1 0];
                
                arrayfun(@(a)... 
                    patch([l(a) l(a) u(a) u(a) ],[0 1 1 0],'yellow','EdgeColor','none', ...
                        'LineWidth',1,'FaceAlpha',.2), ...
                        [1:length(l)])                                
                title('Lower of Bounds');
                subplot(4,1,4); 
                title('Upper of Bounds');
            end
        end
        function res = kdeAnalysis2( obj,  boundType, fwhm, k, showplots )
            % function uses MATLAB Build in KDE
            
            % data: load results array
            gedBoundRes = obj.getGedBoundsRes;
            
            % data: branch by bound type
            ivector = gedBoundRes.ivector.(boundType);

            % results: parameters
            frex = gedBoundRes.frex;   % base frequency range
            upscale = 2; % upsample multipler for KDE (higher more smooth)
            
            % results: preallocate
            kde             = struct();                   % results array
            frex_kde        = linspace(min(frex), ...     % upscale frex 
                                max(frex),length(frex) * upscale);
            clust_kde       = zeros(1, length(frex_kde)); % clust #
            ivector_kde     = zeros(1, length(frex_kde)); % kde results
            impulse_kde     = zeros(1, length(frex_kde)); % kde results
            
            % results: calculate
            impulse         = frex( find(ivector) );      % bound frex
            index           = dsearchn(frex_kde', ...     % match index
                                impulse');                %  to upscale
            impulse_kde(index) = true; 
            frex_kde(index) = impulse;                    % bound vec
            clust_kde(index)= find(impulse);              % clust # vec
            
            pd = fitdist(impulse',...                  % create filt
                'Kernel',...                              % kernal
                'Kernel','Normal',...                     % ~1-1.5 Hz FWHM
                'Width', .4);
            
            ivector_kde = pdf(pd, frex_kde);          % kdensity
            % ivector_kde = mat2gray(ivector_kde);
             figure(1); plot(frex_kde, ivector_kde);
            switch boundType
                case 'lower'               
                    obj.gedBoundsRes.ivectorKde.loclust     = clust_kde; 
                    obj.gedBoundsRes.ivectorKde.lofrex      = frex_kde; 
    
                case 'upper'
                    obj.gedBoundsRes.ivectorKde.upclust     = clust_kde; 
                    obj.gedBoundsRes.ivectorKde.upfrex      = frex_kde; 
            
            end
            obj.gedBoundsRes.ivectorKde.(boundType)        = impulse_kde;
            obj.gedBoundsRes.ivectorKde.filtKde.(boundType)= ivector_kde;
            
            obj.gedBoundsRes.ivectorKde.frex = frex_kde;
            
        end     
        function res = kdeAnalysis( obj,  boundType, fwhm, k, showplots )
            
            if nargin < 2
                fprintf(['\n\nError: Please select a bound type:\n', ...
                    '\n1. range', ...
                    '\n2. edge', ...
                    '\n3. lower', ...
                    '\n4. uppper\n\n'] );
                
                assert(nargin > 2, 'Not enough input arguments');
            end
            
            if nargin < 3
                % full-width half-maximum: the key Gaussian parameter
                fwhm = 4;
                k = 1/10;
                fprintf(['\n\nGenerate KDE Vectors\nUsing defaults FWHM 4 (~1-2Hz)' ,...
                    'and window of 1/10\n\n']);
                showplots = false;
            end
            
            gedBoundRes = obj.getGedBoundsRes;
            frex = gedBoundRes.frex;
            
            switch boundType
                case 'range'
                    ivector = gedBoundRes.ivector.range;                    
                case 'edge'
                    ivector = gedBoundRes.ivector.edge;                    
                case 'lower'
                    ivector = gedBoundRes.ivector.lower;
                case 'upper'
                    ivector = gedBoundRes.ivector.upper;
                    
                otherwise
            end
              
            % FWHM = pd.sigma*2*(2*log(2))^(1/2)
            
            % create and implement Gaussian window
            
            % impulse vector in points
            k = ceil((1/10) * length(ivector));
            sWindow = -k:k;
            
            % create Gaussian window ( 2.355 * sigma = fwhm )
            gauswin = exp( -(4*log(2)*sWindow.^2) / fwhm^2 );
            gausPre = gauswin;
            gauswin = gauswin / sum(gauswin); % normalization
            
            zero_padding_width = ceil(length(frex) * 0.5);
            % zero_padding_width = 50;
            pad_ivector = padarray(ivector,[0 zero_padding_width], 0);
            assert(numel(ivector) == numel(pad_ivector(zero_padding_width+1:length(ivector)+ ...
                zero_padding_width)))
            filtsigG = zeros(1,size(pad_ivector,2));
            
            % implement the weighted running mean filter
            for i=k+1:length(pad_ivector)-k-1
                filtsigG(i) = sum( pad_ivector(i-k:i+k).*gauswin );
            end
            filtsigG = filtsigG(zero_padding_width+1:...
                length(ivector)+ ...
                zero_padding_width);
            
            obj.gedBoundsRes.ivector.filt.(boundType) = filtsigG;
            obj.gedBoundsRes.ivector.fwhm = fwhm;
            
            if showplots == true
                figure(2);
                subplot(3,2,1:2); plot(frex, ivector);
                tstr = sprintf('%s (%s) Original Impulse Vector (0 1) of gedBounds', ...
                    obj.s.subj_basename, obj.s.subj_subfolder);
                title(tstr, 'Interpreter', 'none');
                xlabel('Frequency (Hz)');
                subplot(3,2,3);
                plot(frex(1:length(gausPre)), gausPre); title('Before normalization');
                xlabel('frex(1:length(gauswin)');
                subplot(3,2,4);
                plot(frex(1:length(gauswin)), gauswin); title('After normalization');
                xlabel('frex(1:length(gauswin)');
                subplot(3,2,5:6); plot(frex, filtsigG);
                title(sprintf('Gaussian weighted (FWHM=%d) running mean filter', fwhm));
                xlabel('Frequency (Hz)'); hold on;
                plot(frex, PDF,'r');
                % save the image
                fname = fullfile(obj.pathdb.figures, ...
                    [obj.s.subj_basename '_gedFilt.png']);
                saveas(gcf,fname);
                clf
            end
            
        end
        function fig = plotBandBoundries( obj )
            
            gedBoundsRes = obj.gedBoundsRes;
            
            freqbands = gedBoundsRes.freqbands;
            evecCorMat = gedBoundsRes.evecCorMat;
            frex = gedBoundsRes.frex;
            evecs = gedBoundsRes.evecs;
            evals = gedBoundsRes.evals;
            filterType = gedBoundsRes.filterType;
            maps = gedBoundsRes.maps;
            
            nbchan = num2str(length(maps));
            subTitle = sprintf('%s (%s) %s %s', ...
                obj.s.subj_basename, ...
                obj.s.subj_subfolder, ...
                filterType,...
                nbchan);
            
            %% correlation matrix and band boundaries
            % Cohen
            % complete
            figure(1), clf, colormap(flipud(bone))
            
            imagesc(evecCorMat), hold on
            f2u = round(linspace(1,length(frex),10));
            set(gca,'clim',[0 1],'xtick',f2u,'xticklabel',...
                round(frex(f2u),1),'ytick',f2u,'yticklabel',round(frex(f2u),1))
            axis square, axis xy
            ylabel('Frequency (Hz)');
            
            hold on;
            for i=1:max(freqbands)
                
                tbnds = frex(freqbands==i);
                tbnds = dsearchn(frex',tbnds([1 end])');
                
                for vi = 1 : length(tbnds)
                    box
                    plot(tbnds,[1 1]*tbnds(1),'m','linew',2)
                    plot(tbnds,[1 1]*tbnds(2),'m','linew',2)
                    plot([1 1]*tbnds(1),tbnds,'m','linew',2)
                    plot([1 1]*tbnds(2),tbnds,'m','linew',2)
                end
                
            end
            yyaxis right
            ylim([0 2])
            yticks([0:.1:1])
           % ylabel('Frequency (Hz)'); 
            e = evals(:,1);
            e = e-min(e); 
            e=e./max(e);
            plot(e,'Color',[0, 0, 1, 0.5],'linew',1)
            ylim([0 2]);
            xlabel('Frequency (Hz)');
            ylabel('Normalized Eigenvalue');
            h = colorbar;
            colormap(flipud(bone));
            %set(h,'Ticks',0:.2:1,'fontsize',15) % 'ticklabels',{'1','.8','.6','.4','.2'},
            %set( h, 'YDir', 'reverse' );
            set(get(h,'label'),'string','Eigenvector Correlation (r^2)',...
                'fontsize',10);

            subTitle = sprintf('%s (%s) %s', ...
                obj.s.subj_basename, ...
                obj.s.subj_subfolder, ...
                obj.gedBoundsRes.filterType);
            sgtitle(subTitle, 'Interpreter', 'none');
            
            fname = fullfile(obj.pathdb.figures, ...
                [obj.s.subj_basename '_' nbchan '_' filterType '_bandBounds.png']);
            fname2 = fullfile(obj.pathdb.figures, ...
                [obj.s.subj_basename '_' nbchan '_' filterType '_bandBounds.fig']);
          
            obj.filedb.bounds = fname;
            saveas(gcf,fname);
            saveas(gcf,fname2);

            close gcf
            
        end
        function obj = plotBoundsAndTopo( obj )
            
            obj.plotBandBoundries;
            obj.plotPcaTopography;
            f1 = obj.filedb.topoplot;
            f2 = obj.filedb.bounds;
            obj.combineTwoImages(f1, f2);
            
        end
        function fig = plotPcaTopography(obj)
            gedBoundsRes = obj.gedBoundsRes;
            freqbands = gedBoundsRes.freqbands;
            maps = gedBoundsRes.maps;
            chanlocs = gedBoundsRes.chanlocs;
            frex=gedBoundsRes.frex;
            filterType = gedBoundsRes.filterType;
            
            nbchan = num2str(length(maps));
            subTitle = sprintf('%s (%s) %s %s', ...
                obj.s.subj_basename, ...
                obj.s.subj_subfolder, ...
                filterType,...
                nbchan);
            %% plot the average maps
            figure(2);
            clf;
            colormap jet;
            
            tmpcsv = {};
            
            for i=1:max(freqbands)
                subplot(6,5,i)
                [m,score,latent, ~, explained, ~] = pca(maps(freqbands==i,:));
               % size(mean(maps(freqbands==i,:),1))
                topoplotIndie(m(:,1),...
                    chanlocs,'numcontour',0,...
                    'electrodes','off','plotrad',.6);
                perVar = round((max(latent) / sum(latent)) *100);
                title([ num2str(round(mean(frex(freqbands==i)),2)) 'Hz %var(' num2str(perVar) ')' ])
                
                tmpcsv{i,1} = obj.s.subj_basename;
                tmpcsv{i,2} = i;
                for ci = 1 : length(chanlocs)
                    tmpcsv{i,2+ci} = m(ci,1);
                end
            end
            
            obj.topoRow = tmpcsv;
            
            
            sgtitle(subTitle, 'FontSize',14);
            fname = fullfile(obj.pathdb.figures, ...
                [obj.s.subj_basename '_' nbchan '_' filterType '_topo.png']);
            fname2 = fullfile(obj.pathdb.figures, ...
                [obj.s.subj_basename '_' nbchan '_' filterType '_topo.fig']);
            
            obj.filedb.topoplot = fname;
            saveas(gcf,fname);
            saveas(gcf,fname2);
            close gcf
            
        end
        function fig = plotSpectrogram( obj, bandrange )
            
            if nargin < 2
                bandrange = [2 90];
            else
            end
            
            figure;
            EEG = obj.getCurrentEEG;
            
            figure; pop_spectopo(EEG, 1, [], ...
                'EEG' , 'percent', 15, 'freq', [6 10 22 55], 'freqrange',bandrange,...
                'electrodes','on');
            fname = fullfile(obj.getPermutationDirectory(), ...
                [obj.s.subj_basename '_spectroPo.png']);
            obj.clearLastEEG;
            saveas(gcf,fname);
            close gcf
        end
        function fig = combineTwoImages( obj, f1, f2 )
            [p, fn, ext] = fileparts(f1);
            fn2 = [obj.s.subj_basename '_topoBounds.png'];
            A = imread( f1 );
            B = imread( f2 );
            combImg = imfuse(A, B, 'montage');
            J = imresize(combImg,.75);
            % imshow(J);
            newfilename = fullfile(p, fn2);
            imwrite(J,newfilename);
        end
    end
    methods (Static)
        
        function timew = getTimeWindowFromEEG( EEG )
            timew = [EEG.xmin EEG.xmax];
        end
    end
end

% select data
% average trials
% multiple co-variance matrix for each trial average

% interpret covar matrix
% underlying strucutre
% strength of the diagnonal or variability

% create covariance matrix

% sort channels to make more meaningful

% visualizing quadratic form


