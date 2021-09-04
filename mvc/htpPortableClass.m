classdef htpPortableClass < handle & htpBaseMaster
    % HTPPORTABLECLASS
    %   Light version of high throughput pipeline
    % for sharing analyses with collaborators.
    
    properties
       % sub;
       % htpcfg;
    end
    
    methods
        function obj = htpPortableClass()
            
        end
    end
    
    methods % data import
        function obj  = importDataObjects( obj, csvfile, matfile, datapath )
            % Import: Select valid subjects via CSV
            T = readtable(csvfile, 'Delimiter', ',');
            validIndex = strcmp(T.proc_state(:), 'postcomps');
            
            % Import: Load data objects
            % Exclusion report in CSV
            load(matfile, 'sub');
            arrayfun(@(sub) sub.setCsv( csvfile ), sub, 'uni', 0);
            arrayfun(@(sub) sub.setMat( matfile ), sub, 'uni', 0);
            
            % Import: Prune objects to valid subjects
            % Objects contain preprocessing details / methods
            sub_valid = sub( validIndex );
            obj.sub = sub_valid;
        end
        function res = exportPreprocessingDetails( obj )
            % in: subject info out: subject CSV + details
            % preprocessing results
            varOfInterest = {...
                'subj_subfolder',...
                'epoch_trials', ...
                'proc_xmax_raw',...
                'proc_xmax_epoch',  ...
                'proc_badchans', ...
                'proc_removeComps'};
            
            subjTable = obj.createSubjectTable();
            
            % preallocation
            if length(obj.sub)>1
                tmptable = cell(length(subjTable), ...
                    numel(varOfInterest)); % original -R 03122021
            else
                tmptable = cell(1,numel(varOfInterest));
            end
            
            % get variable of interest
            for i = 1 : size(tmptable,1)
                s = obj.sub(i);
                tmptable{i,1} = i;
                tmptable{i,2} = s.subj_basename;
                for j = 1 : numel(varOfInterest)
                    varname = varOfInterest{j};
                    switch varname
                        case {'proc_badchans', 'proc_removeComps'} % accounts for proc_ipchans in newer builds
                            if isempty(s.proc_ipchans)
                                tmpvalue = numel(s.(varname));
                            else
                                tmpvalue = s.proc_ipchans;
                            end
                        otherwise
                            tmpvalue = s.(varname);
                            
                    end
                    tmptable{i, j+2} = tmpvalue ;
                end
            end
            
            tableHeader = ['id' 'eegid' varOfInterest];
            
            tmptable = cell2table(tmptable);
            tmptable.Properties.VariableNames = tableHeader;
            res = tmptable;
        end
        function res = createSubjectTable( o )
            
            baseHeaderArr = {'id', 'eegid', 'group'};
            baseCellArr = cell(numel(o.sub), ...
                numel(baseHeaderArr));
            
            for sub_i = 1 : numel( o.sub )
                
                s = o.sub(sub_i);
                
                id = sub_i;
                eegid = s.subj_basename;
                group = s.subj_subfolder;
                
                baseCellArr{sub_i, 1} = id;
                baseCellArr{sub_i, 2} = eegid;
                baseCellArr{sub_i, 3} = group;
                
            end
            res = baseCellArr;
        end
        function o = updateBasePaths(o, basePath)
            sub = o.sub;
            arrayfun(@(sub) sub.setBasepath(basePath), sub, 'UniformOutput', false);
            arrayfun(@(sub) sub.updatePaths(basePath), sub, 'UniformOutput', false);
            o.htpcfg.pathdb = o.sub(1).pathdb;
            o.sub = sub;
            
        end
    end
    
end

