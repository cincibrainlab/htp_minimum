classdef publishClass < handle & publishHelperClass
    % PUBLISHCLASS
    %   high level functions for reproducible research results
    %   an instance can be created and saved for a single manuscript

    properties
        am; % analysisMaster object
        paper;
        data;
        tables;
        figures;
        precompute;
        l; % logger class
        status;
        h; % publishClass helper methods
        pcfg;
    end

    methods % start up methods

        function o = publishClass()
            % constructor
            o.paper = struct();
            o.data = struct();
            o.tables = struct();
            o.figures = struct();
            o.status = struct(); % status codes
            o.precompute = struct();

            o.assignUpdateStatus('init');
            % External log utility (Log4m)
            % {'ALL','TRACE','DEBUG','INFO','WARN','ERROR','FATAL','OFF'}
            o.l = log4m.getLogger('log4m_temp.txt');

            h = publishHelperClass;
        end
        
        function index = getRandomSubgroup( o, num_per_group )
            grpInfo = o.am.getGroupInfo;
            index = [num_per_group * size( grpInfo.idx, 1 ),1];
      
                for i = 1 : size( grpInfo.idx, 1 )
                    grpIdx = find(grpInfo.idx(i,:));
                    selectIndexes = grpIdx(ceil(rand(num_per_group,1)* length(grpIdx)));
                    index = [index selectIndexes];
                end

        end
        function index = getRandomTrials( o, num_trials )
           
            
        end

        function [result, o] = assignUpdateStatus(o, code)

            if nargin > 1

                switch code
                    case 'init'
                        o.status.createAnalysisMasterSuccess = false;
                    case 'createAnalysisMasterSuccess'
                        o.status.createAnalysisMasterSuccess = true;
                    otherwise
                        o.l.error('updateStatus', 'Incorrect Status Code');
                end

            else

                result = false;
            end

        end

        function [result, o] = getUpdateStatus(o, code)

            if nargin > 1

                switch code
                    case 'createAnalysisMasterSuccess'
                        result = o.status.createAnalysisMasterSuccess;
                    otherwise
                        o.l.error('updateStatus', 'Incorrect Status Code');
                end

            else
                result = false;
            end

        end

        function o = createAnalysisMaster(o, ~)
            % generate AnalysisMaster object
            if nargin > 1
                refresh = true;
            else
                refresh = false;
            end

            try
                try
                tmppath = o.data.basepath; % check valid basepath
                catch
                    tmppath = o.pcfg.datapath;
                    o.l.info('Warning','o.data.basepath not set, using alternative.');
                end
            catch exception
                o.l.error('createAnalysisMaster', 'No basepath specified.');
                throw(exception);
            end

            if refresh

                try

                    o.am = htpAnalysisMaster;
                    o.am.init_path;
                    o.am.init_environment;
                    o.am.htpcfg.basePath = tmppath;
                    o.l.info('createAnalysisMaster', 'Success.');
                    o.assignUpdateStatus('createAnalysisMasterSuccess');
                catch exception
                    o.l.error('createAnalysisMaster', 'Error creating AnalysisMaster');
                    throw(exception);
                end

            end

            if o.getUpdateStatus('createAnalysisMasterSuccess')
                % get data from csv
                o.am.setStageCSV(o.pcfg.csvfile, o.pcfg.datapath);
                o.am.updateBasePaths(o.pcfg.datapath);
                o.l.info('createAnalysisMaster', sprintf('CSV: %s', o.pcfg.csvfile));

            end

        end

        function o = assignInitialParameters(o, cfg)
            
            newFields = fields(cfg);            
           
            for i = 1 : numel(newFields)
                o.pcfg.(newFields{i}) = cfg.(newFields{i});
            end
            
            currentFields = fields(o.pcfg);
            o.l.info('pcfg fields', sprintf('[%s] ', currentFields{:}));
            
       
%             try
%                 o.l.info('assignInitialParameters', 'Checking Variables');
%                 o.data.basepath = fullfile(cfg.datapath);
%                 o.l.info('Datapath Set As', o.data.basepath);
%                 o.data.csvfile = cfg.csvfile;
%                 o.l.info('Dataset CSV File', cfg.csvfile);
%                 o.paper.longtitle = cfg.longtitle;
%                 o.l.info('Long Title', o.paper.longtitle);
%                 o.paper.shorttitle = cfg.shorttitle;
%                 o.l.info('Short Title', o.paper.shorttitle);
%                 o.precompute = cfg.precompute;
% 
%             catch exception
%                 str = sprintf(...
%                     '\ncfg.datapath, cfg.longtitle, cfg.shorttitle, cfg.csvfile');
%                 o.l.error('Check for missing cfg fields:', str);
%                 throw(exception);
%             end
        end

    end

    methods % general analysis methods

        function result = checkIfPrecomputeAvailable(o, description)

            if nargin > 1
                availablePrecomputations = cfg.precompute;

                switch description
                    otherwise
                        result = false;
                end

            else

            end

        end
        
        function tbl = exportStudyVariables( o )
            varOfInterest = {'subj_subfolder','epoch_trials', ...
                'proc_xmax_raw','proc_xmax_epoch', 'proc_badchans','proc_removeComps'};
            
            o.tables.studyvars = createBaseTable( o.am );
            if length(o.am.sub)>1
                tmptable = cell(length(o.tables.studyvars), ...
                    numel(varOfInterest)); % original -R 03122021
            else
                tmptable = cell(1,numel(varOfInterest));
            end
                        
            % get variable of interest
%             for i = 1 : length(tmptable)
            for i = 1 : size(tmptable,1)
                tmptable{i,1} = i;
                tmptable{i,2} = o.am.sub(i).subj_basename;
                for j = 1 : numel(varOfInterest)
                    varname = varOfInterest{j};
                    switch varname
                        case {'proc_badchans', 'proc_removeComps'} % accounts for proc_ipchans in newer builds
                            if isempty(o.am.sub(i).proc_ipchans)
                                tmpvalue = numel(o.am.sub(i).(varname));
                            else
                                tmpvalue = o.am.sub(i).proc_ipchans;
                            end
                        otherwise
                            tmpvalue = o.am.sub(i).(varname);
                        
                    end
                    tmptable{i, j+2} = tmpvalue ;
                end
            end
            
            tableHeader = ['id' 'eegid' varOfInterest];
            
            tmptable = cell2table(tmptable);
            tmptable.Properties.VariableNames = tableHeader;
            [~,csvbasename,~] = fileparts(o.am.htpcfg.csvfile);
            savefile = o.addPathToFile('analysis',[csvbasename '_QiVars.csv']);

            writetable(tmptable,savefile);
            tbl = tmptable;
            o.l.info('PublishClass:exportStudyVariables', savefile);

        end

        function filename = savePrecomputeData(o, cfg) 

            type = cfg.type;
            filename = cfg.savefile;
            bandTable = cfg.bandTable;
            
            switch type                                  
                case 'power'
                    restlabels = {'freqTable', 'pntsTable', 'rest_abs_power', 'rest_rel_power', 'rest_abs_power', ...
                        'rest_abs_hz', 'rest_rel_hz', 'rest_abs_power_band_average_trials', ...
                        'rest_rel_power_band_average_trials', 'rest_peakloc'};
                case 'relpow'
                    restlabels = {'freqTable', 'pntsTable', 'rest_rel_power', ...
                        'rest_rel_hz', ...
                        'rest_rel_power_band_average_trials', 'rest_peakloc'};
                case 'rest_dbwpli'
                    restlabels = {'rest_dbwpli'};
                    
                case 'networkFeatures'
                    restlabels = {'degree', 'density', 'strength', 'cluster_coef1', ...
                        'cluster_coef2', 'transitivity1', 'transitivity2', 'modularity1', 'modularity2', 'efficiency1', 'efficiency2',...
                        'assortativity1', 'assortativity2', 'betweencentrality1', 'betweencentrality2', 'eigencentrality1', 'eigencentrality2'};
                    try
                        filename = [cfg.networkLabel '_' cfg.savefile];
                        cfg.savefile = filename;
                        %cfg.savefile.sigElectrodesOnly = cfg.savefile_sigElectrodesOnly;
                        disp(['Using significant electrode mode:' filename]);
                    catch
                        %filename = cfg.savefile_sigElectrodesOnly;
                        disp(['Using all electrode mode:' filename]);
                    end
                otherwise
            end
            
            
            [baseCellArr baseHeaderArr]  = createBaseTable( o.am );
            
            switch type
                case 'networkFeatures'
                    % create header and placeholder table
                    headerArray = {'id','eegid'};
                    cellTable = cell(numel(o.am.sub), ...
                        numel(headerArray) + numel(restlabels)*numel(bandTable));
                    counter = numel(headerArray) + 1;
                    
                    for iLabel = 1:numel(restlabels)
                        label = restlabels{iLabel};
                        for iBand = 1 : numel(bandTable)
                            colLabel = [ bandTable{iBand} '_' label ];
                            headerArray{counter} = colLabel;
                            counter = counter + 1;
                        end
                    end
                    tmpvalue = struct();
                    rowNum = 1;
                    for iSubject = 1:numel(o.am.sub)
                        s = o.am.sub(iSubject);
                        cellTable{iSubject, 1} = rowNum;
                        cellTable{iSubject, 2} = s.subj_basename;
                        colNum = 3;
                        for iLabel = 1:numel(restlabels)
                            label = restlabels{iLabel};                            
                            for iBand = 1 : numel(bandTable)
                                tmpFeatureArray = s.(label);                                
                                cellTable{iSubject, colNum} = tmpFeatureArray(iBand);
                                colNum = colNum + 1;
                                rowNum = rowNum + 1;
                            end
                        end
                        
                    end
                    
                    lefttable = cell2table( baseCellArr, 'VariableNames', baseHeaderArr );
                    innertable = cell2table( cellTable, 'VariableNames', headerArray );
                    alltable = outerjoin(lefttable,innertable,'Type','left','MergeKeys',true, 'Keys', {'eegid','eegid'});
                
                    cfg.table = alltable;
                    o.exportTable( cfg );
                    noAdditionalExport = true;
                case 'relpow'
                    tmpvalue = struct();
                    
                    for iSubject = 1:numel(o.am.sub)
                        s = o.am.sub(iSubject);
                        
                        for iLabel = 1:numel(restlabels)
                            label = restlabels{iLabel};
                            tmpvalue(iSubject).(label) = s.(label);
                            [~, fn, ext] = fileparts(filename);
                            s.filename.(['precompute_' type]) = [fn ext];
                        end
                        
                    end
                                        noAdditionalExport = false;

                                       
                otherwise
                    
            end
            if ~noAdditionalExport
                filename =  fullfile(s.pathdb.analysis, filename);
                save(filename, 'tmpvalue');
                disp(['Saved file: ' filename]);
            end
        end

        function savefile = exportTable(o, cfg)

            if nargin < 2
                availableTypes = {'rest_peakloc','networkFeatures'};
                o.l.warn('exportTable, available types', sprintf('%s', availableTypes{:}));
            else

                % check inputs
                try
                    type = cfg.type;                    
                    savefile = o.addPathToFile('analysis', cfg.savefile);
                catch exception
                    o.l.error('exportTable', 'check CFG inputs.');
                    throw(exception);
                end

                switch type
                    
                    case 'networkFeatures'
                        table = cfg.table;
                        writetable( table, savefile );
                        
                    case 'rest_peakloc'

                        % check inputs for this case

                        try

                            selectedElectrodes = cfg.indexScalpChannels;
                            chanLabels = {cfg.chanlocs(:).labels};
                            chanLabels = chanLabels(selectedElectrodes);
                        catch e
                            o.l.error('exportTable:rest_peakloc', 'check CFG inputs.');
                            throw(e);
                        end

                        otherColumnsOfInterestLabels = {'subj_basename', 'subj_subfolder'};
                        otherMeasuresOfInterestLabels = {'Mean_All_Elec', 'Mean_Scalp_Elec'};
                        nSub = numel(o.am.sub);
                        peakloc_row = [];
                        final_row = {};
                        global_count = 1;

                        for iSub = 1:nSub

                            s = o.am.sub(iSub);

                            for iCol = 1:numel(otherColumnsOfInterestLabels)
                                peakloc_row{iSub, iCol} = s.(otherColumnsOfInterestLabels{iCol});
                            end

                            tmprow = s.('rest_peakloc');
                            allelectrodes = {nanmean(tmprow)};

                            selectrodes = {nanmean(tmprow(selectedElectrodes))};

                            final_row(iSub, :) = horzcat(peakloc_row(iSub, :), allelectrodes, selectrodes, num2cell(tmprow));
                            global_count = global_count + 1;
                        end

                        colNames = {};
                        headerLabels = [otherColumnsOfInterestLabels otherMeasuresOfInterestLabels chanLabels];

                        for iHeader = 1:numel(headerLabels)
                            colNames{iHeader} = headerLabels{iHeader};
                        end
                        
                        tbl = cell2table(final_row);
                        try
                            tbl.Properties.VariableNames = colNames;
                        catch e
                            if strcmp('MATLAB:table:VariableNameNotValidIdentifier', e.identifier)
                                tbl.Properties.VariableNames = genvarname(colNames);
                            end
                        end
                        writetable(tbl, savefile);
                        o.l.info('exportTable', sprintf('CSV exporting %s to %s', type, savefile));
                    
                    case 'rest_dbwpli_csv'
                        
                        nBands = size( o.am.sub(1).rest_dbwpli, 3);
                        idxSelChans = cfg.indexScalpChannels;
                                                distarr = o.pcfg.distarr;
                                                
                        % fn in: am.sub(x).rest_dbwpli: chan x chan x band x sub
                        % fn out: chan x chan per band
                        figure;
                        grpinfo = o.am.getGroupInfo;
                            grpidx = grpinfo.idx;
                        counter = 1;
                        for iBand = 1 : nBands
                            conMat = o.getConMat( {o.am.sub(:).rest_dbwpli}, iBand );
                            conMat = conMat( idxSelChans, idxSelChans, :);
                            
                            threshold(iBand) = o.getConnThreshold( conMat, 1 );
                            
                            tband = threshold(iBand); % threshold per band
%                             f = o.plotDistributionThreshold( o.getVectorFromUpperTriMatrix(conMat), threshold(iBand) );
%                             titlestr = sprintf('dbwPLI Band: %s Threshold: %1.3f', o.pcfg.BandLabels2{iBand},   threshold(iBand) );
%                             title(titlestr)
%                             
%                             savefile = fullfile(o.am.htpcfg.pathdb.figs, ['PLI_Threshold_Band_'  o.pcfg.BandLabels2{iBand} '.png']);
%                             o.saveAndTitleFigure(f, savefile);
%                              close(gcf);
                             
                            hubnessChanArr{iBand} = o.calcHubnessFromConnMatrix( conMat, tband );
                            [dist, conn] = o.calcDbwpliByDistance( conMat, tband, distarr(idxSelChans,idxSelChans) );
                            tmphub = zscore(hubnessChanArr{iBand});
                            tmphub = permute(tmphub,[2 1]);
                            
                            grp1hub = mean(tmphub(:, grpidx(1,:)),2);
                            grp2hub = mean(tmphub(:, grpidx(2,:)),2);
                            grp1dist = mean(conn(:,grpidx(1,:)),2);
                            grp2dist = mean(conn(:,grpidx(2,:)),2);
                           % figure; f =  scatter(dist(:,1),grp1dist); hold on; scatter(dist(:,1),grp2dist); title(sprintf('Distance vs. dbwPLI: %s', o.pcfg.BandLabels2{iBand}));
                           
                           clear g
                           
                           g(1,1)=gramm('x',dist(:,1),'y', {grp1dist;grp2dist},  'color', {'FXS';'TDC'});
                           
                           g(1,1).geom_point();
                           
                           g(1,1).set_names('x','Distance (cm)','y','dbwPLI (phase synchrony)');
                           g(1,1).stat_ellipse('type','95percentile','geom','area','patch_opts',{'FaceAlpha',0.1,'LineWidth',2});
                           
                           
                           g(1,1).set_title(sprintf('Distance vs. dbwPLI: %s', o.pcfg.BandLabels2{iBand}));
                           figure('Position',[100 100 800 800]);
                           
                           g.draw();
                           f = gcf;
                           savefile = fullfile(o.am.htpcfg.pathdb.figs, ['PLI_vs_Distance_'  o.pcfg.BandLabels2{iBand} '_grps.png']);
                           o.saveAndTitleFigure(f, savefile);
                           close(gcf);
                           chanlocs = cfg.chanlocs(cfg.indexScalpChannels);
%                             
%                           
%                             subplot(2, nBands,counter); topoplotIndie(grp1hub, chanlocs, 'numcontour', 0); 
%                             title(['FXS-' o.pcfg.BandLabels2{iBand} '-T' num2str(tband)]);   colormap hot; c1= colorbar;
%                             t = gcf; cbarlimits = t.Children(end-1).Limits;
%                             subplot(2,nBands,counter+nBands); topoplotIndie(grp2hub, chanlocs, 'numcontour', 0); 
%                             title(['TDC-' o.pcfg.BandLabels2{iBand} '-T' num2str(tband)]); colormap hot; c2 = colorbar;
%                             t=gcf; t.Children(end-1).Limits = cbarlimits;
%                             c2.Limits = c1.Limits;
                             counter =counter+1;
                            %                             f = o.plotDistributionThreshold( o.getVectorFromUpperTriMatrix(conMat(:,:, grpinfo.idx(1,:)')), threshold(iBand) );
%                             titlestr = sprintf('Grp1: dbwPLI Band: %s Threshold: %1.3f', o.pcfg.BandLabels2{iBand},   threshold(iBand) );
%                             title(titlestr)
%                             
%                             savefile = fullfile(o.am.htpcfg.pathdb.figs, ['PLI_Threshold_Band_'  o.pcfg.BandLabels2{iBand} '_grp1.png']);
%                             o.saveAndTitleFigure(f, savefile);
%                              close(gcf);
%                            
%                             f = o.plotDistributionThreshold( o.getVectorFromUpperTriMatrix(conMat(:,:, grpinfo.idx(2,:))), threshold(iBand) );
%                             titlestr = sprintf('Grp2: dbwPLI Band: %s Threshold: %1.3f', o.pcfg.BandLabels2{iBand},   threshold(iBand) );
%                             title(titlestr)
%                             
%                             savefile = fullfile(o.am.htpcfg.pathdb.figs, ['PLI_Threshold_Band_'  o.pcfg.BandLabels2{iBand} '_grp2.png']);
%                             o.saveAndTitleFigure(f, savefile);
%                             
%                             close(gcf);
                        end
                        
                    
                        
                        
                        
                                                           
                        allpli = nonzeros(allvalues(:,:));
                        threshold = median(allpli) + std(allpli);
                        figure; clf
                        subplot(311), hold on;
                        
                        histogram(allpli,100);
                        plot([1 1]*threshold, get(gca, 'ylim'), 'r--', 'linew', 3);
                        xlabel('dbwPLI', ylabel('Count'));
                        legend({'Distribution'; 'Threshold'})
                        
                        pliallThresh = allpli > threshold;
                        size(allpli(pliallThresh))
                        
                        grpinfo = o.am.getGroupInfo;
                        
                        grp1 = nonzeros(allvalues(grpinfo.idx(1,:)', :));
                        grp2 = nonzeros(allvalues(grpinfo.idx(2,:)', :));
                        
                        figure; hist([grp1 grp2], 100); hold on;
                        hist(grp2, 100);
                        axis square;
                        figure; imagesc(mean(squarevalues(:,:,grpinfo.idx(1,:)'),3))
                        
                    otherwise
                        o.l.warn('exportTable', 'no export type identified.');

                end

            end

        end
    
        function fileWithPath = addPathToFile(o, type, filename)
            
            try
                fileWithPath = fullfile(o.am.htpcfg.pathdb.(type), filename);
            catch e
                o.l.error('fileWithPath', 'Error with assignment. Check if object loaded.');
                throw(e)
            end
            
        end
        
    end

    methods % spectral methods

        function o = precomputeSpectralPower(o, cfg)
            % inputs required:
            % savefile = precompute savefile
            %
            if nargin < 2
                o.l.error('precomputeSpectralPower', 'Not enough inputs');
            end

            try

                freqBands = cfg.freqBands;
                savefile = fullfile(o.am.htpcfg.pathdb.analysis, cfg.savefile);

            catch exception
                o.l.error('precomputeSpectralPower', 'Missing required cfg input.');
                throw(exception);
            end

            for i = 1:numel(o.am.sub)

                s = o.am.sub(i);

                s.loadDataset('postcomps');

                s.setFreqTable(freqBands); % new 11/22
                s.getPntsTable;
                s.generateStoreRoom;
                s.correlation_global_coupling;
                s.peakdetect_v1;
                s.bandAverageTrials;
                s.subj_trials = s.EEG.trials;
                s.unloadDataset;

            end
            
            cfg.type = 'relpow';
            o.savePrecomputeData( cfg );

        end

        function o = statRelativePower(o, cfg)
            suptitle = @(x) title(x);

            % check inputs
            try
                savefile = o.addPathToFile('analysis', cfg.savefile);
                select_chan = cfg.indexScalpChannels;
                sChanlocs = cfg.sChanlocs;
                nbband = cfg.nBands;
                bandTable = cfg.BandLabels;
                nbgroups = cfg.nbgroups;
                nbchan = cfg.nChans;
                sub = o.am.sub;
                tmp = o.am.sub(1).rest_rel_power_band_average_trials;
                chan_hood = cfg.chanhood;
                individualPlots = cfg.turnOnIndividualPlots;
                groupPlots = cfg.turnOnGroupPlots;
                exportCsv = cfg.turnOnCsvExport;

                clear tmp;
            catch e
                o.l.error('statRelativePower', 'Not enough inputs');
                throw(e);
            end

            newChanlocs = sChanlocs;

            cum_power_rel_band = NaN(nbchan, nbband, length(sub), nbgroups);

            g1 = 0;
            g2 = 0;
            clear list1, clear gender1
            clear list2, clear gender2

            for i = 1:length(sub)
                s = sub(i);

                if ~isempty(s.rest_rel_power_band_average_trials)% change back to rel
                    %         if strcmp(s.rest_subj_gender, 'F')
                    pxx_rel_band = s.rest_rel_power_band_average_trials(:, select_chan)';

                    if strcmp(s.subj_subfolder, o.pcfg.groupInfo.grouplabels{1})%&& strcmp(s.rest_subj_gender, 'M')
                        gid = 1;
                        cum_power_rel_band(:, :, i, gid) = pxx_rel_band;
                        g1 = g1 + 1;
                        list1(g1, 1) = str2double(strtok(s.subj_basename, {'D', '_'}));
                        %    gender1{g1,1} = s.rest_subj_gender;
                    elseif strcmp(s.subj_subfolder, o.pcfg.groupInfo.grouplabels{2})%&& strcmp(s.rest_subj_gender, 'F')
                        gid = 2;
                        cum_power_rel_band(:, :, i, gid) = pxx_rel_band;
                        g2 = g2 + 1;
                        list2(g2, 1) = str2double(strtok(s.subj_basename, {'D', '_'}));
                        %   gender2{g2,1} = s.rest_subj_gender;
                    end

                    %         end
                end

            end

            cum_power_rel_band_g1 = cum_power_rel_band(:, :, ~isnan(cum_power_rel_band(...
                randi(nbchan), randi(nbband), :, 1)), 1); % FXS
            cum_power_rel_band_g2 = cum_power_rel_band(:, :, ~isnan(cum_power_rel_band(...
                randi(nbchan), randi(nbband), :, 2)), 2); % TDC
            mean_power_rel_band_g1 = mean(cum_power_rel_band_g1, 3);
            mean_power_rel_band_g2 = mean(cum_power_rel_band_g2, 3);
            
            [rel_pval1, rel_t_orig1, ~, rel_est_alpha1] = clust_perm2_rui3(...
                cum_power_rel_band_g1, cum_power_rel_band_g2, chan_hood, 5000, .05, 1, .05, 2, [], 1);
            rel_index1 = zeros(nbchan, nbband);
            rel_chan1 = NaN(nbchan, nbband);
            sig_gt_chan_rel = cell(1, nbband);

            [rel_pval2, rel_t_orig2, ~, rel_est_alpha2] = clust_perm2_rui3(...
                cum_power_rel_band_g1, cum_power_rel_band_g2, chan_hood, 5000, .05, -1, .05, 2, [], 1);
            rel_index2 = zeros(nbchan, nbband);
            rel_chan2 = NaN(nbchan, nbband);
            sig_ls_chan_rel = cell(1, nbband);

            for m = 1:nbband
                rel_index1(:, m) = logical(rel_pval1(:, m) < rel_est_alpha1(m));
                rel_chan1(:, m) = select_chan' .* rel_index1(:, m);
                sig_gt_chan_rel{m} = nonzeros(rel_chan1(:, m));
                rel_index2(:, m) = logical(rel_pval2(:, m) < rel_est_alpha2(m));
                rel_chan2(:, m) = select_chan' .* rel_index2(:, m);
                sig_ls_chan_rel{m} = nonzeros(rel_chan2(:, m));
            end

            if individualPlots
                % Individual figures
                savepath = o.am.htpcfg.pathdb.figs;

                tp = @(x) topoplot(x, newChanlocs, 'headrad', 'rim', 'electrodes', 'on', 'gridscale', 300, 'hcolor', 'k');
                cb = @() colorbar('SouthOutside');

                bandno = [1:7];

                for i = 1:length(bandno)
                    handle = figure;
                    m = bandno(i);

                    lower_bound = min(min(mean_power_rel_band_g1(:, m)), min(mean_power_rel_band_g2(:, m)));
                    upper_bound = round(max(max(mean_power_rel_band_g1(:, m)), max(mean_power_rel_band_g2(:, m))), 1);

                    subplot(2, 1, 1)

                    [h, grid_or_val, plotrad_or_grid, xmesh, ymesh] = tp(mean_power_rel_band_g1(:, m));
                    caxis([0, upper_bound])
                    %h=colorbar('SouthOutside','Ticks', 0:upper_bound/3 :upper_bound);

                    subplot(2, 1, 2)

                    tp(mean_power_rel_band_g2(:, m));
                    caxis([0, upper_bound])

                    %h=colorbar('SouthOutside', 'Ticks', 0:upper_bound/3 :upper_bound);

                    ha = findobj(handle, 'type', 'axes');

                    for i = 1:length(ha)
                        ha(i).XLimMode = 'auto';
                        ha(i).YLimMode = 'auto';
                    end

                    saveas(gcf, fullfile(savepath, ['B', num2str(m), '.png']));
                    %fig2svg(fullfile(savepath, ['B' num2str(m) '.svg']), gcf);

                    close(handle);
                end

            end

            % Full Figure with comparisons
            if groupPlots
                figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            end

            for m = 1:nbband
                lower_bound = min(min(mean_power_rel_band_g1(:, m)), min(mean_power_rel_band_g2(:, m)));
                upper_bound = max(max(mean_power_rel_band_g1(:, m)), max(mean_power_rel_band_g2(:, m)));

                if groupPlots
                    subplot(4, nbband, m)
                    topoplot(mean_power_rel_band_g1(:, m), newChanlocs, 'headrad', 0)
                    upper_bound = round(upper_bound, 1);
                    caxis([0, upper_bound])

                    h = colorbar('northoutside', 'Ticks', [0, upper_bound]);

                    % fig1=figure;
                    % left=100; bottom=100 ; width=20 ; height=500;
                    % pos=[left bottom width height];
                    % axis off
                    % caxis([0 upper_bound])
                    % colorbar;
                    % set(fig1,'OuterPosition',pos)
                    % savepath = am.htpcfg.pathdb.figs;
                    % saveas(fig1, fullfile(savepath, ['colorbar' num2str(m) '.png']));
                    % close(fig1);
                    %

                    set(h, 'Position', [.088 * (m - 1) + .145, .715, .04, .013]); %9 bands
                    %     set(h, 'Position', [.115*(m-1)+.155 .715 .04 .013]);  % 7 bands

                    title(bandTable(m), 'FontSize', 16)

                    subplot(4, nbband, m + nbband)
                    topoplot(mean_power_rel_band_g2(:, m), newChanlocs, 'headrad', 0)
                    caxis([0, upper_bound])
                end

            end

            dim = [.05, .66, .2, .18]; % x y w h
            str = {['FXS (', num2str(g1), ')']};
            annotation(gcf, 'textbox', dim, 'EdgeColor', 'none', ...
                'String', str, 'FitBoxToText', 'on', 'FontSize', 12);
            dim = [.05, .45, .2, .18]; % x y w h
            str = {['TDC(', num2str(g2), ')']};
            annotation(gcf, 'textbox', dim, 'EdgeColor', 'none', ...
                'String', str, 'FitBoxToText', 'on', 'FontSize', 12);

            for m = 1:nbband
                rel_index1(:, m) = logical(rel_pval1(:, m) < rel_est_alpha1(m));
                rel_t_orig_sig1(:, m) = rel_t_orig1(:, m) .* rel_index1(:, m);

                if groupPlots
                    subplot(4, nbband, m + 2 * nbband);
                    topoplot(rel_t_orig_sig1(:, m), newChanlocs, 'headrad', 0);
                    caxis([-5, 5])
                end

                if m == nbband
                    h = colorbar('SouthOutside');
                    set(h, 'Position', [.088 * (m - 1) + .145, .29, .04, .013]);
                    %         set(h, 'Position', [.115*(m-1)+.155 .29 .04 .013]); % 7 bands
                end

                rel_index2(:, m) = logical(rel_pval2(:, m) < rel_est_alpha2(m));
                rel_t_orig_sig2(:, m) = rel_t_orig2(:, m) .* rel_index2(:, m);

                if groupPlots
                    subplot(4, nbband, m + 3 * nbband);
                    topoplot(rel_t_orig_sig2(:, m), newChanlocs, 'headrad', 0);
                    caxis([-5, 5])
                end

            end

            if groupPlots
                dim = [.04, .235, .2, .18]; % x y w h
                str = {'FXS > TDC'};
                annotation(gcf, 'textbox', dim, 'EdgeColor', 'none', ...
                    'String', str, 'FitBoxToText', 'on', 'FontSize', 12);

                dim = [.04, .025, .2, .18]; % x y w h
                str = {'FXS < TDC'};
                annotation(gcf, 'textbox', dim, 'EdgeColor', 'none', ...
                    'String', str, 'FitBoxToText', 'on', 'FontSize', 12);

                spt = suptitle('P1-dipole FXS vs TDC relative power');
                spt.FontSize = 16;

                [~, savefile_fig, ~] = fileparts(savefile);
                savefile_fig = o.addPathToFile('figs', [savefile_fig '.png']);
                saveas(gcf, savefile_fig);
            end

            if exportCsv
                %% Publication Figures
                sigelec1 = rel_t_orig_sig1;
                sigelec2 = rel_t_orig_sig2;
                sigelec1_log = sigelec1(:, :) ~= 0;
                sigelec2_log = sigelec2(:, :) ~= 0;

                % arrays with signifigant channels
                % Publication Data Output for Correlations

                % col 1: subj_basename
                c1 = {sub.subj_basename}';
                c2 = {sub.subj_subfolder}';

                avgpow = squeeze(mean(cum_power_rel_band, 1));

                catarr = zeros(1, size(avgpow, 1));

                for i = 1:size(avgpow, 3)

                    tmparr = avgpow(:, :, i)';
                    tmparr = squeeze(tmparr(~isnan(tmparr(:, 1)), :));

                    catarr = vertcat(catarr, tmparr);

                    tmparr = [];
                end

                catarr(1, :) = [];

                t = table(c1, c2);
                t2 = array2table(catarr);
                final_table = [t, t2];

                avgpow = squeeze(mean(cum_power_rel_band, 1));

                all_sigelc = horzcat(sigelec1_log, sigelec2_log);
                all_subs = cat(3, cum_power_rel_band_g1, cum_power_rel_band_g2);

                for i = 1:size(all_sigelc, 2)

                    eleidx(:, i) = logical(all_sigelc(:, i));

                end

                powmat = zeros(size(all_subs(:, :, 1)))
                powtable_sig = zeros(size(all_subs, 3), size(powmat, 2) * 2);

                for i = 1:size(all_subs, 3)

                    powmat(:, :) = all_subs(:, :, i);

                    powrow = zeros(1, size(powmat, 2) * 2);

                    for j = 1:size(powmat, 2)

                        powcol(:, 1) = powmat(:, j);
                        powidx1(:, 1) = eleidx(:, j)';
                        powidx2(:, 1) = eleidx(:, j + 7)';

                        powidx1_n = length(find(powidx1));
                        powidx2_n = length(find(powidx2));
                        %disp(powidx1_n)

                        if powidx1_n ~= 0, p1 = mean(powcol(logical(powidx1)));
                            else, p1 = NaN;
                        end

                        if powidx2_n ~= 0, p2 = mean(powcol(logical(powidx2)));
                            else, p2 = NaN;
                        end

                        powrow(j) = p1;
                        powrow(j + 7) = p2;

                        fprintf('B: %d/%d E1: %d(%1.3f) E2: %d(%1.3f)\n', j, j + 7, powidx1_n, p1, powidx2_n, p2)

                    end

                    % final power table with just significant electrodes
                    powtable_sig(i, :) = powrow;

                end

                powtable_sig_t = array2table(powtable_sig);
                alltable = horzcat(final_table, powtable_sig_t);

                bandTable = ["delta", "theta", "alpha1", "alpha2", "beta", "gamma1", "gamma2"];

                varNames = bandTable;
                bandTablecell = arrayfun(@(x)char(varNames(x)), 1:numel(varNames), 'uni', false);
                varNames = {'eegid', 'group'};

                for i = 1:length(bandTablecell)

                    fxsp{i} = [bandTablecell{i}, '_fxs'];
                    tdcp{i} = [bandTablecell{i}, '_tdc'];

                end

                alltable.Properties.VariableNames = [varNames, bandTablecell, fxsp, tdcp];

                writetable(alltable, savefile);
            end

        end

    end

    methods % connectivity methods

        function precomputeDebiasedPhaseLagIndex(o, cfg)
            % check inputs
            % check inputs

            try
                usePrecomputedFile = cfg.usePrecomputedFile;
                precomputedSaveFile = o.addPathToFile('analysis', cfg.precomputedSaveFile);

                if ~usePrecomputedFile
                    savefile = o.addPathToFile('analysis', cfg.savefile.dbwpli);
                end

            catch e
                o.l.error('statRelativePower', 'Not enough inputs');
                throw(e);
            end
            
            if strcmp(cfg.datatype,'signal')
               datatype = 'signal'; 
            else
               datatype =  'postcomps';
            end
            
            if ~usePrecomputedFile

                % calculate measure for each subject
                N = numel(o.am.sub);
                sub = o.am.sub;

                for i = 1:N
                    s = sub(i);
                    tic;
                    fprintf('Processing: %d/%d %s', i, N, s.subj_basename);
                    s.loadDataset(datatype);
                    if strcmpi(datatype,'signal')
                        s.study_type = 'rest';
                        s.epochData;
                    end
                    s.ram_conn_dwpli;
                    s.unloadDataset;
                    toc;
                    sub(i) = s;
                end

                % save results to variable
                rest_dwpli = {};

                for i = 1:N
                    s = sub(i);
                    rest_dwpli{i} = s.rest_dbwpli;

                end

                save(fullfile(o.am.htpcfg.pathdb.analysis, ['rest_dwpli_60s.mat']), 'rest_dwpli');

            else

                o.loadPreComputedData('rest_dbwpli', precomputedSaveFile);

            end

        end

        function o = statGroupDebiasedPhaseLagIndex(o, cfg)
            % check inputs
            
            try
                usePrecomputedFile = cfg.usePrecomputedFile;
                if usePrecomputedFile
                    precomputedSaveFile = o.addPathToFile('analysis', cfg.precomputedSaveFile);
                end
                select_chan = cfg.indexScalpChannels;
                chanlocs = cfg.chanlocs;
                distArrFile = cfg.distArrFile;
                load(distArrFile, 'distarr'); 
                chanlocsMniFile = cfg.chanlocsMniFile;
                nbband = cfg.nBands;
                bandTable = cfg.BandLabels2;
                nbgroups = cfg.nbgroups;
                chan_hood = cfg.chanhood;
                savefile_cluster = o.addPathToFile('analysis', cfg.savefile);

            catch e
                o.l.error('statGroupDebiasedPhaseLagIndex', 'Not enough inputs');
                throw(e);
            end
            % dbwpli = ch x ch x band
            sub = o.am.sub;

            dist_threshold = 10;

            pliArr = zeros(length(chanlocs), size(sub(1).rest_dbwpli, 3), numel(sub));

            for sub_i = 1:numel(sub)

                s = sub(sub_i);

                for pow_i = 1:size(s.rest_dbwpli, 3)

                    mat = s.rest_dbwpli(:, :, pow_i);
                    plimat = triu(mat) + triu(mat, 1)';

                    for chan_i = 1:length(plimat)
                        chpli = plimat(:, chan_i);
                        d = distarr(:, chan_i);

                        idx = d > dist_threshold;

                        ch.mean = mean(chpli(idx));
                        ch.sd = std(chpli(idx));
                        ch.median = median(chpli(idx));
                        ch.max = max(chpli(idx));
                        ch.min = min(chpli(idx));

                        pliArr(chan_i, pow_i, sub_i) = ch.mean;

                    end

                end

            end

            grpInfo = o.am.getGroupInfo;
            bandpool = [1:nbband];
            
            %bandname = {'Delta', 'Theta', 'Alpha1', 'Alpha2', 'Beta', 'Gamma1', 'Gamma2'};
            bandname = bandTable;
            
            idxMF = {o.am.sub(:).subj_gender};
            
            grp1idx = grpInfo.idx(1,:);     % group 1 index
            grp2idx = grpInfo.idx(2,:);     % group 2 index
            
            grp1s = idxMF( grp1idx );
            grp2s = idxMF( grp2idx );
            
            % grp1s_idx = strcmp('M', grp1s);
            % grp2s_idx = strcmp('M', grp2s);
            
            grp1_n = length(find(grp1idx));
            grp2_n = length(find(grp2idx));
            
            % grp1idx = logical([grp1s_idx zeros(1, grp2_n)]);
            % grp2idx = logical([zeros(1, grp1_n) grp2s_idx]);
            
            
            
            grp1allbands = squeeze(pliArr(:,:, grp1idx)); % individual subjects
            grp2allbands = squeeze(pliArr(:,:, grp2idx));
            
            size(grp1allbands)
            
            % starting with chan / sub / band
            % target is chan / band / sub
            
            grp1 = grp1allbands; %permute(grp1, [1 3 2]);
            grp2 = grp2allbands; %permute(grp2, [1 3 2]);
            
            % scalp channels only
            grp1 = grp1(select_chan, :, :);
            grp2 = grp2(select_chan, :, :);
            
            % Cluster based permutation testing
            if ~usePrecomputedFile
                
            [pval1, t_orig1, ~, est_alpha1] = clust_perm2_rui3(grp1, grp2, ...
                chan_hood(select_chan, select_chan), 5000, 0.05, 1, .05, 2, [], 1); % FXS>TDC
            
            [pval2, t_orig2, ~, est_alpha2] = clust_perm2_rui3(grp1, grp2, ...
                chan_hood(select_chan, select_chan), 5000, 0.05, -1, .05, 2, [], 1); % FXS<TDC
            
            precompute.dbwpli.pval1       = pval1;
            precompute.dbwpli.t_orig1     = t_orig1;
            precompute.dbwpli.est_alpha1  = est_alpha1;
            precompute.dbwpli.pval2       = pval2;
            precompute.dbwpli.t_orig2     = t_orig2;
            precompute.dbwpli.est_alpha2  = est_alpha2;
                        
            save( savefile_cluster, 'precompute');
            
            else
                load( precomputedSaveFile, 'precompute')
                
                pval1       = precompute.dbwpli.pval1;
                t_orig1     = precompute.dbwpli.t_orig1;
                est_alpha1  = precompute.dbwpli.est_alpha1;
                pval2       = precompute.dbwpli.pval2;
                t_orig2     = precompute.dbwpli.t_orig2;
                est_alpha2  = precompute.dbwpli.est_alpha2;
                
            end
            
            size(t_orig1)
            colorvec = {};
            
            % plotting
            grpnames = {'FXS', 'TDC'};
            
            grp1all = pliArr(:,:, grp1idx); % individual subjects
            grp2all = pliArr(:,:, grp2idx);
            
            
            grpnames1 = repmat({sprintf('FXS (n=%d)', grp1_n) }, 1, 108);
            grpnames2 = repmat({sprintf('TDC (n=%d)', grp2_n) }, 1, 108);
            
            grpnameall = [grpnames1 grpnames2];
     
            grp1mean = mean(grp1all, 3);   % mean subjects
            grp2mean = mean(grp2all, 3);
            
            bardat = [grp1mean(select_chan,:)' grp2mean(select_chan,:)'];
            %bardat = zscore(bardat);
            
          
%             clear g;
%             figure;
%             for iBand = 1 : size( bardat, 1 )
%                 g(1,iBand) = gramm('x', grpnameall, 'y', bardat(iBand,:), 'color', grpnameall);
%                 g(2,iBand) = g(1,iBand).copy;
%                 g(1,iBand).stat_summary('geom',{'bar','black_errorbar'});
%                 %g(1,1).axe_property('YLim',[0.035 0.045]);
%                 g(1,iBand).set_title(bandname{iBand});
%                 g(2,iBand).geom_jitter('width',0.4,'height',0);
%                 g(1,iBand).no_legend();
%                 g(2,iBand).no_legend();
%                 % g(1,2).set_title('geom_jitter()');
%                 %figure('Position',[100 100 800 550]);
%                 g.set_names('x','Group','y','dbwPLI','color','Group');
%                 g.set_title(sprintf('Long Distance (10 cm) Phase Synchrony'));
%                 
%             end
%             g.draw();
%             saveas(gcf, o.addPathToFile('figs', ['dbwpli_long_electrode_all_bands.png']));  
%             close(gcf);
            
            g1 = struct();
            g2 = struct();
            g1.nameidx = '';
            g2.nameidx = '';
            g1.bandidx = '';
            g2.bandidx = '';
            
            g1pli = [];
            g2pli = [];
            
            
            
            % Electrode Average Subplots
            for band_i = 1 : size(t_orig1,2)
                
                bandname{band_i};
                g1pli = [g1pli mean(grp1(:,band_i,:),3)'];
                g2pli = [g2pli mean(grp2(:,band_i,:),3)'];
                
                g1.length = length(mean(grp1(:,band_i,:),3)');
                g2.length = length(mean(grp1(:,band_i,:),3)');
                
                g1.nameidx = [g1.nameidx repmat({'FXS'}, 1, g1.length)];
                g2.nameidx = [g2.nameidx repmat({'TDC'}, 1, g2.length)];
                
                g1.bandidx = [g1.bandidx repmat({bandname{band_i}}, 1, g1.length)];
                g2.bandidx = [g2.bandidx repmat({bandname{band_i}}, 1, g2.length)];
                
                %colorvec = [colorvec bandname{band_i} bandname{band_i}];
                
            end
            
            figure('Position',[100 100 1600 800]);
            clear g;
            g(1,1) = gramm('x', [g1.nameidx g2.nameidx], 'y', [g1pli g2pli], 'color', [g1.nameidx g2.nameidx]); %, 'subset', strcmp('Delta', [g1.bandidx g2.bandidx]) & strcmp('Beta', [g1.bandidx g2.bandidx])  );
            g(2,1) = g(1).copy;
            g(1,1).geom_jitter('width',0.4,'height',0);
            g(1,1).facet_grid([],[g1.bandidx g2.bandidx], 'scale','independent');
            g(2,1).facet_grid([],[g1.bandidx g2.bandidx], 'scale','independent');
            g(2,1).stat_boxplot();
            g(1,1).no_legend();
            g(2,1).no_legend();
            g(1,1).set_order_options('column',0);
            g(2,1).set_order_options('column',0);
            
            g(1,1).set_names('x','','y','dbwPLI','column','');
            g(2,1).set_names('x','Group','y','dbwPLI','column','');
            g.set_title('Mean Long Range Phase/Phase Synchrony per Electrode by Frequency Band')
            g.draw();
            
            saveas(gcf, o.addPathToFile('figs', ['dbwpli_long_electrode_all_bands.png']));

            
        end
        
        function o = statGroupDbwpliNBS( o, cfg )
            
            % check and assign inputs
            
            try
                indexScalpChannels = cfg.indexScalpChannels;
                nBands = cfg.nBands;
                nbchan = cfg.nChans;
                bandTable = cfg.BandLabels2;
                alpha = cfg.alpha;
                exe = cfg.surficeExecutable;
                networkFeaturesCsv = cfg.networkFeaturesCsv;
                contrastType = cfg.contrastType;  % tail of test
                grpInfo = o.am.getGroupInfo;
                surficeMNI = o.am.htpcfg.chanNow.net_surficeMni;
                surfice_on = cfg.surficeOutputOn;
                mnielec = readtable(surficeMNI);
                restlabels = cfg.networkFeatures;

            catch e
                o.l.error('statGroupDbwpliNBS', 'Error in CFG assignments.');
                throw(e);
            end
            
            % initalize variables                   
            cum_dwpli = []; n = []; conn = []; p = [];
            
            % create skeleton 'base' subject table
            sig_elec_tbl_header = {'gid', 'group', 'eegid'};
            sig_elec_tbl = {};
            
            [baseCellArr baseHeaderArr]  = createBaseTable( o.am );
            lefttable = cell2table( baseCellArr, 'VariableNames', baseHeaderArr );

            % load( fullfile(o.am.htpcfg.pathdb.analysis, cfg.precomputedSaveFile ) );
            
                        
            for iBand = 1 : nBands
                
                contrastTypeArray = [1 2];      % tails for significance testing
                
                for iContrast = 1 : numel(contrastTypeArray)
                    
                    contrastType = contrastTypeArray( iContrast );                    
                    
                    [cum_dwpli, n, conn, p, GLM] = htp_dwpli_nbs(o.am, iBand, contrastType, alpha);
                    
                     if ~isempty(conn)                         
                        
                        for iConn = 1 : numel(conn)
                            
                            c = conn{iConn};
                            edgefile= o.addPathToFile('analysis', ['NBS_PLI_' bandTable{iBand} '_' num2str(iConn) '_.txt']);
                            edge = full(c);
                            
                            % create CSV of measures of this node
                            
                            if networkFeaturesCsv == true
                                % FUNCTION: SAVE SIGNIFICANT NETWORK PROPERTIES                                
                                    o.am.sub = o.clearNetworkFeaturesFromCsv( o.am.sub, restlabels );
                                    o.am.sub = o.computePliNetworkFeatures( o.am.sub, indexScalpChannels, edge ); % prune select_chan to edge only

                                    cfg.savefile_sigElectrodesOnly = ...
                                        o.generatePliSigElectrodeCsvFilename( o, edgefile );
                                
                                    cfg.bandTable = { bandTable{ iBand } }; %#ok<CCAT1>
    
                                    % save network measures                                
                                    o.savePrecomputeData( cfg );
        
                            end
                            
                            
                            if surfice_on == true
                                
                                % generate figure
                                cfg.edge = edge;
                                cfg.savefileEdgeFile = edgefile;                                
                                o.generateSurficeVisualization(o, cfg)
                                
                            end
                            
                            % generate CSV
                            % get all subjects
                            cfg.iBand = iBand;
                            cfg.contrastType = contrastType;
                            cfg.edgefile = edgefile;
                            
                            resultColumnWithKey = o.createPliSubCsvColumn(o, cfg, o.am.sub, edge );
                            
                            lefttable = outerjoin(lefttable,resultColumnWithKey,'Type','left','MergeKeys', true, 'Keys', {'eegid','eegid'});                            
                            
                        end
                        
                        
                        
                    end
                    
                end
                
  
            end
            
            sig_elec_tbl = lefttable;
            saveTable = sig_elec_tbl;
           % saveTable.Properties.VariableNames = sig_elec_tbl_header;
            [~,csvfile, ~] = fileparts( cfg.savefile );
            csvfile = fullfile(o.am.htpcfg.pathdb.analysis, [csvfile '_.csv']);
            writetable( saveTable, csvfile );
            
            
%         c = table2cell(mnielec);
%         c = c(pcfg.indexScalpChannels',:);
%         x = ''; y=''; z = ''; in = ''; co='';
%         for i = 1 : length( c )
%             
%             X = c{i,1};
%             Y = c{i,2};
%             Z = c{i,3};
%             
%             IN = 1.1;
%             CO = 3;
%             
%             x = sprintf('%s,%2.1f',x,X);
%             y = sprintf('%s,%2.1f',y,Y);
%             z = sprintf('%s,%2.1f',z,Z);
%             
%             in = sprintf('%s,%2.1f',in,IN);
%             co =  sprintf('%s,%2.1f',co,CO);
%         end
%         
%         bracketstr = @(x) sprintf('[%s]', strip(x,'Left',','));
%         
%         tmpX = bracketstr(x);
%         tmpY = bracketstr(y);
%         tmpZ = bracketstr(z);
%         tmpin = bracketstr(in);
%         tmpco = bracketstr(co);
%         
%         cmd.nodecreate = sprintf('NODECREATE('''', %s, %s, %s, %s, %s);', ...
%             tmpX, tmpY, tmpZ, tmpin, tmpco);



        end
        
    end
    
    methods % graph theory methods
        
    end
    
    methods (Static)
        
        function [dist, conn] = calcDbwpliByDistance( conMat, tband, distArr )
            
            for i = 1 : size(conMat,3)
                
               v(:,i) = publishClass.tri2vec(conMat(:,:,i));
               d(:,i) = publishClass.tri2vec(distArr);
              % figure; scatter(d,v)
            end
            dist = d;
            conn = v;
            
        end
        
        function v = tri2vec( mat )
            At = mat(:,:,1).';
                m = (1:size(At)).' > (1:size(At,2));
                v = At(m);
        end
        function conMat = getConMat( conMatFull, band )

            size(conMatFull); % 1 x 138
            size(conMatFull{1});
            
            conMat = NaN( size(conMatFull{1},1),  size(conMatFull{1},2), numel(conMatFull) );
            % 128 x 128 x 138
            
            for i = 1:numel(conMatFull)
              
                tmpdat = conMatFull{i};
                conMat(:,:,i) = tmpdat(:,:, band);
                
            end
            
             
        end
        
        function hubnessChanArr = calcHubnessFromConnMatrix( conMat, threshold )
            for i = 1 : size(conMat,3)
                tmpConMat(:,:,i) = (triu(conMat(:,:,i)) + triu(conMat(:,:,i),1)');
                thresConMat(:,:,i) = tmpConMat(:,:,i) > threshold(end);
                hubnessChanArr(i,:) = sum(thresConMat(:,:,i))/(size(conMat,1)-1);
            end   
        end
        
        function threshold = getConnThreshold( conMatAllSubByBand, nSd )
            % conMatAllSubByBand 128 x 128 x 138 (subjects)
            % nSd = number of standard deviations + median
            
            allvalues = zeros(numel(nonzeros(triu(conMatAllSubByBand(:,:,1),1))), size( conMatAllSubByBand, 3));
            
            for iSub = 1 : size( conMatAllSubByBand, 3)
            
                allvalues(:, iSub) = nonzeros(triu( conMatAllSubByBand(:,:,iSub), 1 ));
            
            end
            
            threshold = median(nonzeros(allvalues)) + (nSd * std(nonzeros(allvalues)));
        end
        
        function vec = getVectorFromUpperTriMatrix( mat )
            
            for i3 = 1 : size( mat, 3)
                
                vec(:, i3) = nonzeros(triu( mat(:,:,i3), 1 ));
                
            end
            vec = nonzeros(vec);
        end
        
        function saveAndTitleFigure( fig, savefile )
            % title(fig, sprintf('%s', titlestr));
            saveas(fig, savefile, 'png');
            
        end
        function fig = plotDistributionThreshold( mat, threshold)
            
            fig = figure; histogram(nonzeros(mat), 50); hold on;
            plot([1 1]*threshold, get(gca, 'ylim'), 'r--', 'linew', 3);
            xlabel('dbwPLI'), ylabel('Count');
            legend({'Distribution'; 'Threshold'});
           
            return;
        end
        
        function mat = data2array( arr )
            
        end
        
    end

end
