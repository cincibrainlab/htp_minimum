classdef publishHelperClass < handle
    % PUBLISHHELPERCLASS Smaller methods and helper functions
    
    properties
        
    end
    
    methods
        function o = publishHelperClass()
            
        end
    end
    
    methods % deal with precomputed files
        function assignDataPath( o, datapath )
        
            o.pcfg.datapath  = fullfile(datapath); 
            o.data.basepath = o.pcfg.datapath;
        end
        
        function assignSubjectCsv( o, csvfile_subject )
        
            o.pcfg.csvfile  = csvfile_subject; 
            o.data.csvfile = o.pcfg.csvfile;
        end
        
        function assignChanLocsMat( o, chanlocsmat )
            o.pcfg.chanlocsMat = chanlocsmat;
        end
        
        function chanlocs = getChanLocs(o)
            chanlocs = load(o.pcfg.chanlocsMat, 'chanlocs');
            chanlocs = chanlocs.chanlocs;
            o.pcfg.chanlocs = chanlocs;
            o.pcfg.ChanlocsN = numel(chanlocs);

        end
        
        function chanhood = getChanhood( o )
            chanhood = load(o.pcfg.chanhoodMat);
            chanhood = chanhood.chan_hood;
            o.pcfg.chanhood = chanhood;
            
        end
        
        function chanlocs = getScalpChanlocs( o, indexScalpChannels )
            
            if nargin < 2
                indexScalpChannels = o.pcfg.indexScalpChannels;
            end
            
            chanlocs = o.getChanLocs;
            o.pcfg.sChanlocs = chanlocs( indexScalpChannels );
            chanlocs = o.pcfg.sChanlocs;
            o.pcfg.sChanlocsN = numel(chanlocs);
        end
        
        function assignScalpChannels( o, idx )
            o.pcfg.indexScalpChannels = idx;
        end
        function chanhood = getScalpChanhood( o, indexScalpChannels )
            if nargin < 2
                indexScalpChannels = o.pcfg.indexScalpChannels;
            end
            chanhood = o.getChanhood;
            o.pcfg.sChanhood = chanhood( indexScalpChannels, indexScalpChannels );
            
        end
        
        function assignBandDefs( o, bandDefs, bandDefsExtended )
            if nargin > 2
            
                o.pcfg.BandLabels = bandDefs(:,1);
                o.pcfg.nBands = numel(o.pcfg.BandLabels);
                o.pcfg.freqBands = cell2mat(bandDefs(:,2:3));
                
                bandDefs = [bandDefs; bandDefsExtended]; 
                o.pcfg.BandLabelsExtended = bandDefs(:,1);
                o.pcfg.freqBandsExtended = cell2mat(bandDefs(:,2:3));
                o.pcfg.nBandsExtended = numel(o.pcfg.BandLabelsExtended);
                
            else
                
                o.pcfg.BandLabels = bandDefs(:,1);
                o.pcfg.freqBands = cell2mat(bandDefs(:,2:3));
                o.pcfg.freqBands = cell2mat(bandDefs(:,2:3));
            end
        end
        
        function assignGroupInfo( o )
            o.pcfg.groupInfo = o.am.getGroupInfo;
            o.pcfg.nGroups = o.pcfg.groupInfo.nbgroups;
        end
        
        function pathdb = getPathDb( o )
            pathdb = o.am.htpcfg.pathdb;
        end
        
        function assignGroupLabels( o, cellarr )
            if numel(cellarr) ~= o.pcfg.nGroups
                o.l.error('Error', ' Not enough labels for number of groups.');
            else
                o.pcfg.groupInfo.grouplabels2 = cellarr;
            end
        end
        
            
            
                
          
        
        function chanLabels = getChanLabels( o, indexScalpChannels )
            chanlocs = o.getChanLocs;
            if nargin < 2                
                chanLabels = {chanlocs.labels};
                o.pcfg.ChanLabels = chanLabels;
                indexScalpChannels = o.pcfg.indexScalpChannels;
                chanlocs = chanlocs(indexScalpChannels);
                chanLabels = {chanlocs.labels};
                o.pcfg.sChanLabels = chanLabels;
            else
                chanlocs = chanlocs(indexScalpChannels);
                chanLabels = {chanlocs.labels};
                o.pcfg.sChanLabels = chanLabels;
            end
            
        end

        
        function assignChanMniTxt( o, chanMniTxt )
            o.pcfg.chanMniTxt = chanMniTxt;
        end
        
        function assignChanhoodMat( o, chanhoodMat )
            o.pcfg.chanhoodMat = chanhoodMat;
        end           
        
        function assignChanDistanceMat( o, chandistMat )
            o.pcfg.chandistMat = chandistMat;
        end
        
        
        function loadPreComputedData( o, cmd, filename )
           % load precomputed data into structures
           % switch commmand
           % filename of precomputed file
           
           switch cmd
               
               case 'clear_relpow'
                   %tmpvar = load(filename);
                   %relpow = tmpvar.tmpvalue;
                   fieldnames = {'freqTable','pntsTable','rest_rel_power','rest_rel_hz',...
                       'rest_rel_power_band_average_trials','rest_peakloc', 'rest_abs_power', ...
                       'rest_abs_hz', 'rest_abs_power_band_average_trials'};
                   
                   for iSub = 1 : numel( o.am.sub )
                       for iField = 1 : numel(fieldnames)
                           o.am.sub(iSub).(fieldnames{iField}) = ...
                               [];
                       end
                   end
                   
               case 'relpow'
                   % full file and path to precomputed relpower in the format:
                   %     struct 1 x number of subjects
                   %     freqTable
                   %     pntsTable
                   %     rest_rel_power
                   %     rest_rel_hz
                   %     rest_rel_power_band_average_trials
                   %     rest_peakloc
                   tmpvar = load(filename);
                   relpow = tmpvar.tmpvalue;
                   fieldnames = {'freqTable','pntsTable','rest_rel_power','rest_rel_hz',...
                       'rest_rel_power_band_average_trials','rest_peakloc'};
                   
                   for iSub = 1 : numel( o.am.sub )
                       for iField = 1 : numel(fieldnames)
                           o.am.sub(iSub).(fieldnames{iField}) = ...
                               relpow.(fieldnames{iField});
                       end
                   end
                   
               case 'rest_dbwpli'
                   tmp = load(filename);
                   loadrest = tmp.rest_dwpli;
                   
                   for i = 1:numel(o.am.sub)
                       o.am.sub(i).rest_dbwpli = zeros(size(loadrest{1, i}));
                       o.am.sub(i).rest_dbwpli = loadrest{1, i};
                   end
                   
               otherwise
               
           end
           
        end
        
        function signal = hideNotchFilter( o, signal, allhz, index )
            specificIndex = dsearchn(allhz', index');
            signal(:, specificIndex(1):specificIndex(2)) = nan;
            
        end
        
        function out = changeLogicalGroupArrToGroupLabels(o, logical_group_array, labeltrue, labelfalse)
            N = numel( logical_group_array );
            
            out = cell(1, N);
            
            for i = 1 : numel(logical_group_array)
               
                if logical_group_array(i)
                    out(i) = {labeltrue};
                else
                    out(i) = {labelfalse};
                end
                
            end
            
        end
        
        function out = createBandLabelArrayFromRanges(o, hz, freqbands, bandlabels )
            
            hz_n = numel(hz);            
            [band_n, ~] = size(freqbands)
            
            % pre-allocate array
            bandCats = cell(1, hz_n);
            bandCats(1:end) = deal({'Undefined'});
            
            for bandi = 1 : band_n
                
                starthz = freqbands(bandi,1);
                endhz = freqbands(bandi,2);
                newLabelIndex = dsearchn(hz', [starthz:endhz]');
                
                bandCats(newLabelIndex(1):newLabelIndex(end)) = deal(bandlabels(bandi));
                
            end
            
            out = bandCats;
            
        end
        
        function f = createChanPowHeatmap( o, chanpow )
            
            grpInfo = o.am.getGroupInfo;
            
            f = figure('position',[800 100 800 1200]);
            
            clf
            mincolor = min([chanpow(:)]);
            maxcolor = max([chanpow(:)]);
            
            grp1idx = grpInfo.idx(1,:)';
            
            subplot(4,1,1)
            heatmap(chanpow(grp1idx,:), o.pcfg.sChanLabels, ...
                1:sum(grp1idx),[],  'TickAngle', 45, 'ShowAllTicks', true, ...
                'Colorbar', true, 'MinColorValue', mincolor, 'MaxColorValue', maxcolor)
            grp2idx = grpInfo.idx(2,:)';
    
            title('WT (n=10) Power per Channel');
            
            subplot(4,1,2)
            heatmap(chanpow(grp2idx,:), o.pcfg.sChanLabels, ...
                1:sum(grp2idx),[],  'TickAngle', 45, 'ShowAllTicks', true, ...
                'Colorbar', true, 'MinColorValue', mincolor, 'MaxColorValue', maxcolor)
            
            title('KO (n=11) Power per Channel');
            
            subplot(4,1,3)
            plot(1:numel(o.pcfg.sChanLabels), [mean(chanpow(grp1idx,:)); mean(chanpow(grp2idx,:))])
            xlabel('Frequency (Hz)');
            ylabel('Rel. Pow');
            legend({'ET','KO'},'Location','eastoutside');
            title('Power by Group per Channel');
            
            chanpow = zscore(chanpow);
            subplot(4,1,4)
            grp1avg = mean(chanpow(grp1idx,:),2);
            x = ones(1, length(grp1avg));
            plot( grp1avg, x, 'b*', 'MarkerSize', 5, 'LineWidth', 2);
            grid off;
            hold on;
            grp2avg = mean(chanpow(grp2idx,:),2);
            
            x = 1.1 * ones(1, length(grp2avg));
            plot(grp2avg, x, 'r*', 'MarkerSize', 5, 'LineWidth', 2);
            grid off;
            hold on;
            % Set up axes.
            ylim([0.5, 1.5]);
            %xlim([mincolor, maxcolor]);
            ax = gca;
            ax.YTick = [1, 1.2];
            ax.YTickLabels = {'WT','KO'};
            xlabel('Power');
            grid on;
            title('Normalized Mean Power Per Subject'); 

        end
        
        
        
    end
    
    methods (Static)
        function [newFilename, fieldname] = generateFileName( type, comment )
            if nargin < 1
               disp('Check inputs.'); 
            end
            
            if nargin <2
                comment = '';
            end
            
            switch type
                case 'relpow'
                    str = 'relpow';
                    ext = 'mat';
                case 'abspow'
                    str = 'abspow';
                    ext = 'mat';
                case 'peakloc'
                    str = 'peakloc';
                    ext = 'csv';
                case 'relpowCsv'
                    str = 'relpow';
                    ext = 'csv';
                case 'dbwpli'
                    str = 'dbwpli';
                    ext = 'mat';
                case 'dbwpliCluster'
                    str = 'dbwpliCluster';
                    ext = 'mat';
                case 'rest_dbwpli_csv'
                    str = 'rest_dbwpli_csv';
                    ext = 'csv';
                case 'networkFeatures'
                    str = 'rest_networkFeatures_csv';
                    ext = 'csv';
                otherwise
            end
            
            r = publishHelperClass.generateRandomNumber();
                
            newFilename = fullfile(['P' r '_' str '_' comment '.' ext]);
            fieldname = [str ext];
        end
        
        function singlerowarray = multirow2singlerow( multirowArr )
            
            tmparray = multirowArr';
            singlerowarray = tmparray(:)';
            
        end
        
        function hz = hz2labels( hz, hzIdx, no_subjects )
            tmphz = hz( hzIdx );
            hz = repmat(tmphz,1 , no_subjects);
            
        end
        function hzGroupLabels = hz2groups( groupLabels, hzIdx)
             hzGroupLabels = repmat(groupLabels, size(hzIdx,2), 1 );
             hzGroupLabels = hzGroupLabels(:)';
            % repmat({'KO', 'WT'}, 1, 4);
        end
        
        function sub = clearNetworkFeaturesFromCsv( sub, restlabels )
            % clear old values
            for iSub = 1 : numel(sub)
                s = sub(iSub);
                for iFeature = 1 : numel(restlabels)
                    s.(restlabels{iFeature}) = [];
                end
                sub(iSub) = s;
            end            
        end
        
        function filename = generatePliSigElectrodeCsvFilename( o, basename )
            
                   [a,b,c] = fileparts( basename );
                    filename = fullfile( [b '.csv'] );

        end
        
        function generateSurficeVisualization(o, cfg)
            
            edgefile = cfg.savefileEdgeFile;
            edge = cfg.edge;
            exe = cfg.surficeExecutable;
            
            plit = array2table(round(edge,2));
            writetable( plit, edgefile, 'Delimiter', '\t','WriteVariableNames', false);
            [a,b,~] = fileparts( edgefile );
            movefile(edgefile,fullfile(a,[b 'FXS.edge']));
            edgefile2 =  fullfile(a,[b 'FXS.edge']);
            basename = fullfile(b);
            
            
            %figname = o.addPathToFile('figs', [basename '.png']);
            node{1} = sprintf('BEGIN RESETDEFAULTS;');
            node{2} = sprintf('MESHLOAD(''mni152_2009''); NODECREATE('''', [69.0,60.8,43.5,32.0,16.3,-0.6,-14.8,49.7,31.0,20.9,0.2,-16.8,-27.9,21.1,0.6,0.5,-20.0,-31.6,-42.5,-17.6,-29.7,-42.7,-49.7,-48.4,-59.7,-62.8,-61.1,-53.6,-42.4,-23.4,-68.0,-74.5,-75.9,-71.9,-62.2,-45.3,-83.0,-83.2,-77.8,-67.7,-84.6,-82.2,-76.6,-77.7,-73.7,-66.6,-48.9,-27.3,-1.6,-76.2,-69.1,-59.7,-48.2,-28.4,-1.9,-63.6,-63.0,-54.1,-41.9,-22.5,-45.8,-41.7,-32.1,-17.4,-2.2,-14.8,-2.2,13.0,18.1,24.3,23.9,20.8,10.9,27.9,37.5,43.9,45.4,42.8,38.2,50.3,55.8,63.3,65.3,43.5,60.5,65.9,71.1,74.3,61.7,74.7,75.5,81.0,76.9,60.8,40.6,12.8,83.5,83.3,71.6,53.0,26.6,83.8,76.6,61.7,42.5,75.0,63.7,50.1], [42.8,56.1,66.9,62.3,48.2,28.1,3.5,71.4,80.6,79.1,72.2,48.2,21.9,85.9,86.1,82.3,79.1,62.4,40.9,86.6,80.8,67.3,50.4,72.4,57.0,42.8,29.3,14.9,-3.5,-21.5,44.2,28.8,15.4,0.6,-11.7,-31.1,-9.7,-17.2,-24.8,-41.1,-40.6,-42.8,-49.8,-64.4,-68.5,-65.8,-59.9,-52.3,-38.1,-63.7,-83.3,-88.6,-83.6,-76.3,-63.5,-73.3,-87.8,-100.8,-101.9,-97.0,-97.8,-107.3,-113.9,-109.2,-102.7,-117.9,-118.7,-109.5,-97.6,-76.9,-53.1,-22.1,-118.5,-114.4,-102.6,-84.5,-61.2,-32.4,-108.4,-102.3,-90.2,-67.6,-43.0,-100.3,-90.0,-85.2,-70.9,-52.0,-75.3,-66.1,-66.7,-45.3,-26.9,-13.2,-4.5,3.1,-43.1,-19.2,-0.9,14.0,21.3,-11.7,14.2,28.7,40.6,27.3,42.2,50.0], [-31.4,0.4,29.5,50.7,72.9,86.1,94.2,-19.4,6.9,29.9,50.8,72.0,82.7,-11.9,10.7,27.9,29.2,49.5,61.3,-12.2,6.3,28.5,44.0,-20.0,-0.7,26.9,46.8,67.7,84.4,97.2,-32.2,-5.5,24.2,50.9,70.5,88.0,-10.0,23.2,50.7,67.6,-5.3,27.1,47.9,0.1,32.2,57.3,81.0,94.6,102.2,-31.1,4.1,36.8,64.0,83.0,95.2,-60.6,-28.3,6.4,38.8,61.4,-54.9,-26.1,7.5,39.6,55.7,-23.1,11.0,39.7,61.6,83.3,95.6,98.1,-22.8,7.9,38.8,64.1,82.1,89.7,-25.4,7.3,37.5,58.6,69.3,-54.6,-27.2,5.2,33.7,49.8,-60.0,-29.8,1.5,29.4,53.5,72.7,85.9,94.3,-3.8,25.8,53.4,69.9,83.4,-8.2,26.6,49.4,63.2,-4.2,28.8,45.7], [1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1], [3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);');
            node{3} = sprintf('EDGELOAD(''%s'')', edgefile2);
            node{4} = sprintf('viewaxial(true);');
            node{5} = sprintf('savebmp(''%s'')',  o.addPathToFile('figs', [basename 'axial.png']));
            
            node{6} = sprintf('viewsagittal(true);');
            node{7} = sprintf('savebmp(''%s'')',  o.addPathToFile('figs', [basename 'sag.png']));
            node{8} = sprintf('quit');
            node{9} = sprintf('END.');
            
            
            tmpscriptname = fullfile(a, [b 'script.txt']);
            fileID = fopen(tmpscriptname,'w');
            for iNode = 1 : numel(node)
                
                fprintf(fileID, '%s\n', node{iNode});
            end
            fclose(fileID);
            
            cmd = sprintf('%s -S %s', exe, tmpscriptname);
            system(cmd);
            
            
        end
        
        function resultColumn = createPliSubCsvColumn(o, cfg, sub, edge )
           
            iBand = cfg.iBand;
            indexScalpChannels = cfg.indexScalpChannels;
            contrastType = cfg.contrastType;
            edgefile = cfg.edgefile;
                        
            [~,col_basename,~] = fileparts(edgefile);
            sig_elec_col_header = [col_basename 'mean_' num2str(contrastType)];
            
            all_dwpli = {sub(:).rest_dbwpli};
            sig_elec_col = {};
            
            for isub2 = 1 : numel(all_dwpli)
                key = sub(isub2).subj_basename;
                current_dwpli = all_dwpli{1, isub2};
                current_dwpli = current_dwpli(indexScalpChannels, indexScalpChannels, iBand);
                sig_edges = current_dwpli( logical(edge) );
                
                sig_elec_col{isub2,1} =  mean(sig_edges);
                sig_elec_col{isub2,2} =  key;
                % sig_elec_tbl{isub2, totalColumns + iConn} = ;
                
            end
            
            resultColumn = cell2table( sig_elec_col, 'VariableNames', [{sig_elec_col_header} 'eegid']);
        end
        
        function sub = computePliNetworkFeatures( sub, scalpchannels, sigEdges )
            
            edge = sigEdges;
            select_chan = scalpchannels;
            nbchan = length(select_chan);
            for i = 1 : length( sub  )
                
                tmp_s = sub(i);
                
                for band = 1 : 1 % delta
                    
                    mat = tmp_s.rest_dbwpli(select_chan, select_chan,band);
                    % remove all other electrodes but
                    % significant network
                    mat( ~logical( edge) ) = nan;
                    
                    % binary
                    thresh_mean = nanmean(mat(:));
                    thresh_std = nanstd(mat(:));
                    bin_mat = mat>thresh_mean-thresh_std & mat<thresh_mean+thresh_std;
                    un_bin_mat = bin_mat + triu(bin_mat)';
                    tmp_s.degree(band) = mean(degrees_und(un_bin_mat)); % degree
                    tmp_s.density(band) = mean(density_und(un_bin_mat)); % density
                    tmp_s.cluster_coef1(band) = mean(clustering_coef_bu(un_bin_mat));
                    tmp_s.transitivity1(band) = mean(transitivity_bu(un_bin_mat));
                    tmp_s.modularity1(band) = mean(modularity_und(un_bin_mat));
                    tmp_s.efficiency1(band) = mean(efficiency_bin(un_bin_mat));
                    tmp_s.betweencentrality1(band) = mean(betweenness_bin(un_bin_mat));
                    tmp_s.eigencentrality1(band) = mean(eigenvector_centrality_und(un_bin_mat));
                    tmp_s.assortativity1(band) = mean(assortativity_bin(un_bin_mat,0));
                    % weighted
                    wu_mat = triu(mat) + triu(mat)';
                    wu_mat(1:nbchan+1:end) = zeros(1,nbchan);
                    % check this with rui, convert to
                    % 0s not nans
                    wu_mat(isnan(wu_mat)) = 0;
                    
                    tmp_s.strength(band) = mean(strengths_und(wu_mat));
                    tmp_s.cluster_coef2(band) = mean(clustering_coef_wu(wu_mat));
                    tmp_s.transitivity2(band) = mean(transitivity_wu(wu_mat));
                    tmp_s.modularity2(band) = mean(modularity_und(wu_mat));
                    tmp_s.efficiency2(band) = mean(efficiency_wei(wu_mat,2));
                    tmp_s.betweencentrality2(band) = mean(betweenness_wei(wu_mat));
                    tmp_s.eigencentrality2(band) = mean(eigenvector_centrality_und(wu_mat));
                    tmp_s.assortativity2(band) = mean(assortativity_wei(wu_mat,0));
                    
                end
                
                sub( i ) = tmp_s;
                
            end
            
        end
               
        function r = generateRandomNumber()
            % generate random number
            a = 10000;
            b = 99999;
            r = round((b - a) .* rand(1, 1) + a);
            r = sprintf('%s', num2str(r));
        end
    end
end

