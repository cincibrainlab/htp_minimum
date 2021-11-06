classdef htpDirectImportClass < handle
    properties
        am   % htpPreprocessMaster
    end

    methods
        function o = htpDirectImportClass()
            o.am = htpPreprocessMaster;
            o.am.firstRun;
        end

        function items = populateChanDropDown(o)
            items = {o.am.xml_elec.net_displayname};
        end
        function items = populateSelectdatafiletypeDropDown(o)
            items = {o.am.xml_elec.net_displayname};
        end

        function o = configHandler( o, key, value )
            switch key
                case 'basePath'
                    o.am.setBasePath(value{1});
                    o.am.study_title = value{2};
                case 'chanInfo' 
                    o.am.htpcfg.chanNow = o.am.htpcfg.chaninfo(value);
                case 'filterInfo'
                    o.am.htpcfg.chanNow.net_filter = o.am.htpcfg.chaninfo(value).net_filter;
            end
        end
        function o = preprocess_direct( o )
            % stage assigment for CSV
            stage_last = 'raw';
            stage_next = 'postcomps';

            % get filelist and subfolders
            o.am.import_filelist;
            fnlist     = o.am.fnlist;
            flength    = length(fnlist);
            o.am.htpcfg.xmax_arr = [];
            o.am.assignConfig2Subjects;

            arrayfun( @(s) s.setElectrodeSystem( o.am.htpcfg.chanNow ), o.am.sub, 'uni', 0);
            arrayfun( @(s) s.setUser( 'direct' ), o.am.sub, 'uni', 0);
            arrayfun( @(s) s.changeStudyTitle( o.am.study_title ), o.am.sub, 'uni', 0);

            totalsubs = length( o.am.sub );
            sub = o.am.sub;

            o.am.setOptNow(repmat(1,[1 25]));
            opt = o.am.formatOptions;

            for i = 1 : totalsubs
                s = sub(i);
                s.EEG = pop_loadset(s.filename.raw, fullfile(s.pathdb.raw, s.subj_subfolder));
                s.EEG = eeg_checkset( s.EEG );
                s.proc_state = 'postcomps';
                s.storeDataset(s.EEG, s.pathdb.postcomps, s.subj_subfolder, s.filename.postcomps);
                s.unloadDataset;
                s.outputRow( stage_next );
                sub(i) = s;
            end

                o.am.sub = sub;

            try
                o.am.createResultsCsv(o.am.sub, 'postcomps','Default');
            catch
            end

        end

    end

end