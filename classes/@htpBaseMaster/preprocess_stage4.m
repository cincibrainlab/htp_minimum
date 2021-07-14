%% Preprocess Stage 4
%
%%% Usage
%   
% obj = preprocess_stage4( obj )
%
%%% Parameters
%
% * INPUTS: obj
%
% * OUTPUTS: obj
%
% The optional input if the function is not self-invoked obj is the 
% htpPreprocessMaster object.  The output is the updated 
% htpPreprocessMaster object with information updated regarding the 
% status and details of stage 4 completion.
%
%%% Copyright and Contact Information
%
% Copyright (C) 2019  Cincinnati Children's (Pedapati Lab)
%
% This file is part of High Throughput Pipeline (HTP)
% 
% See https://bitbucket.org/eped1745/htp_stable/src/master/
% 
% Contact: ernest.pedapati@cchmc.org

%% preprocess_stage4
% Post-ica artifact detection is performed, automatic if the user selects 
% the stage 4 auto button, or manually if the regular stage 4 button is 
% selected, to finish the preprocessing of the group of subjects.
%
% If the user is not on the right stage, it is determined if the subject has
% successfully completed stage 4 prior to the current moment.  If that is
% the case the subject will be saved to the appropriate directory and
% counted as completed.  If the subject is not on the right stage and has
% not previously completed stage 4, it will be marked as an error and the
% file will be skipped.
%
% The data rank, set by the user in the gui settings, determines how many
% initial components are plotted to the user through the selectcomps.m
% function.  The user will have multiple windows of time series, components,
% and plots (power activity, colormaps of , etc.to guide their decision for 
% artifact removal.  The various figures are placed next to the central 
% messagebox that the user can then reject specific components by entering 
% the number of the component they wish to reject.  Once all determined 
% artifacts are removed, then the user can save their work for the subject 
% into the postica directory and continue on to the next subject.  The user
% can also save their current progress and come back to the stage later to 
% finish the in progress subjects. 
%
% The steps taken for stage 4 depend upon the type of data being preprocessed.  
% If the data is Resting data, then the stage will proceed as described
% above.  If there are events in the dataset, then there are various paths
% taken depending on the type of event after the execution of the above 
% component selection process. 
%
% If Event datasets,the subject will undergo stage A of din separation and 
% be prompted to utilize the options dialog application to define the event 
% epochs in the dataset.  Then, the data will be cleaned and saved with 
% different paths taken for stage A epochs and every other stage din with 
% stage A having extra parsing steps if the level b attribute of the event 
% configuration for each subject is set before cleaning and saving.  
%
% However, if the stage is B then level B of din separation will be 
% performed. The other staged dins will proceed with cleaning of the epochs 
% and saving for the subject without first parsing like for stage A. 
function obj = preprocess_stage4( obj )

stage_last = 'postica';
stage_next = 'postcomps';

[mc, mm, mw] = obj.tools_log;       % output log
opt     = obj.formatOptions;        % get GUI options


dipfitcalc = opt.dipfitcalc;

disp_welcome(mm);


obj.sub = obj.loadSub( stage_last );  % retrieve subs from spreadsheet


obj.composeComponentCsv(obj.sub, str2double(opt.pca_rank));


arrayfun(@( s ) s.setopt( opt ), obj.sub, 'uni', 0);  % assign current options to each sub

objStageStatus              = obj.htpcfg.objStageStatus;
objStageStatus_completed    = obj.htpcfg.objStageStatus_completed;

% initalize loops
prev_files = 0; skip_files = 0; errorchk = 0; uniquei = 1; flg = 0;

% 

if strcmp(obj.htpcfg.optnow.Stage2_CleanMode, 'FullAuto') || strcmp(obj.htpcfg.optnow.Stage2_CleanMode, 'ASR')
    
    sub = obj.sub;
    
    for i = 1 : length(sub)
        s = sub(i);
        if ismember(i, objStageStatus)
            s.loadDataset( 'postica' );
            
            % automatically select components
            % algo: out of first 50 components
            % selects cardiac, eye, line, and muscle components
            % that have > threshold confidence (i.e. 0.75 out of 1).
            % excludes file if > 10% of total # of total components
            % for generating auto drafts for testing only            
            s.removeComps_auto( 0.75, false);     
            
            % remove components
            s.compRemove;
            
            % Generic DIPFIT calculations
            if strcmpi(dipfitcalc,'On')
                s.icview;
                s.calc_dipoles;
            else
                s.set_dipfitcalc( 0 );
            end
            
            % save cleaned dataset into the postcomps directory
            
            obj.tool_saveSubject( s, stage_next);
            
            % unload data & decrease memory footprint
            s.unloadDataset;
            
            s.outputRow('postcomps');
            
        end
        
        sub(i) = s;
    end
    
    obj.sub = sub;
    
else
    % manual mode
    while uniquei <= length(obj.sub) && flg ~=2
        
        s = obj.sub(uniquei);   % allows for parallel processing and temporary changes
        
        if ismember(uniquei, objStageStatus) % processes only files in correct stage
            
            s.setUser( obj.user );
            s.changeStudyTitle( obj.study_title );
            
            s.subj_id = [uniquei length(obj.sub)];
            
            switch opt.epoch_type
                
                case 'Event'
                    
                    s.loadDataset( 'postica' );
                    s = tool_selectComps( obj, s);
                    mc(sprintf('Components Selected: %s\n', mat2str(s.proc_removeComps)));
                    %s = s.compRemove;
                    
                    s.setStudyType( opt.epoch_type );
                    if ~s.exclude_switch
                        % Level A of DIN Separation
                        s.eventcfg.levelb = 0;
                        s.eventcfg.fixsettings = 0;
                        
                        if strcmp(opt.epoch_type, 'Event') == 1
                            uihandle = optionsDialogApp( s, obj );
                            uiwait(uihandle.UIFigure);
                        end
                        obj.tool_saveSubject( s, stage_next);
                        s = obj.tool_createEpochs( s );
                        obj.htpcfg.eventcfg = s.eventcfg;
                        
                        % break up file into "big" cuts
                        % select event classifiers
                        % locate epochs that have classifiers within the big
                        % epoch
                        % make array of epochs
                        
                        switch s.eventcfg.stage
                            case 'A'
                                
                                if s.eventcfg.levelb == 1
                                    
                                    for j = 1 : length(s.eventcfg.event_fn)
                                        
                                        s.loadDataset_Event( j );
                                        
                                        s = dinparser( s );
                                        
                                        if strcmp(s.eventcfg.stage, 'C')
                                            % s.eventcfg.details(3:end,:) = [];
                                            for z = 2 : length(  s.eventcfg.event_fn' )
                                                s.loadDataset_Event( z );
                                                
                                                s = obj.tool_cleanEpochs( s );
                                                obj.tool_saveSubject( s, stage_next);
                                    
                                                s.unloadDataset;
                                                s.outputRow('postcomps');
                                            end
                                        end
%                                         s = obj.tool_cleanEpochs( s );
%                                     
%                                         obj.tool_saveSubject( s, stage_next);
%                                     
%                                     
%                                         s.unloadDataset;
%                                         s.outputRow('postcomps');
                                    end
                                    
%                                     s = obj.tool_cleanEpochs( s );
%                                     
%                                     obj.tool_saveSubject( s, stage_next);
%                                     
%                                     
%                                     s.unloadDataset;
%                                     
%                                     s.outputRow('postcomps');
                                    
                                else
                                    for j = 1 : length(s.eventcfg.event_fn)
                                        s.loadDataset_Event( j );
                                        s = obj.tool_cleanEpochs( s );
                                        obj.tool_saveSubject( s, stage_next);
                                        
                                        
                                        s.unloadDataset;
                                        s.outputRow('postcomps');
                                        
                                    end
                                    
                                end
                                
                            case 'B'
                                if s.eventcfg.levelb == 1
                                    s.loadDataset_Event(1);
                                    
                                    % Level B of DIN Separation
                                    if strcmp(opt.epoch_type, 'Event') == 1
                                        uihandle = optionsDialogApp( s, obj );
                                        uiwait(uihandle.UIFigure);
                                    end
                                    obj.tool_saveSubject( s, stage_last);
                                    
                                    s.eventcfg.fixsettings = 1; % set next sub files identical
                                    obj.htpcfg.eventcfg = s.eventcfg;
                                    
                                    
                                    s = obj.tool_createEpochs( s );
                                end
                        end
                        
                        
                    else
                        s.unloadDataset;
                        
                        s.outputRow('postcomps');
                    end
                    
                case 'Rest'
                    
                    s.loadDataset( 'postica' );
                    
                    s = tool_selectComps( obj, s);
                    
                    % save cleaned dataset into the postcomps directory
                    
                    obj.tool_saveSubject( s, stage_next);
                    
                    % unload data & decrease memory footprint
                    s.unloadDataset;
                    
                    s.outputRow('postcomps');
                    
                otherwise
            end
            
            [flg, errorchk] = obj.redoOrContinue( s );
            
        else
            
            if ismember(uniquei, objStageStatus_completed)
                s.unloadDataset;
                s.proc_state = 'postcomps';
                s.outputRow(stage_next);
                prev_files = prev_files + 1;
                flg = 0;
            else
                s.unloadDataset;
                s.outputRow('error');
                skip_files = skip_files + 1;
            end
            
        end
        
        % Code to verify or redo subject
        switch flg
            case 0
                obj.sub(uniquei) = s;  % apply changes to persistant object
                uniquei = uniquei + 1;   % increment counter
            case 1
                % does not increment counter
                obj.msgout( 'Redo subject selected.');
                s = [];
            case 2 % added 9/4 save work EP
                
                mm('Save work in progress...');
                break;
        end
        
    end
end

stagecnt = @( stage_name ) sum(strcmpi( stage_name, {obj.sub(:).proc_state} ));
stagedisp = @( stage_name, n ) sprintf('Stage 4 %s: %d', stage_name, n);

mm( stagedisp( 'Total Valid:', stagecnt( stage_last ) + stagecnt( stage_next ) )  );
mm( stagedisp( stage_last, stagecnt( stage_last ) ) );
mm( stagedisp( stage_next, stagecnt( stage_next ) ) );
[n, idx] = obj.countExcluded;
mm( stagedisp( 'Excluded:', n ));

if ~isempty( idx )
    
    for ei = 1 : length( idx )
        
        s = obj.sub(idx(ei));
        mm( sprintf('  ID: %d\t%s (%d): %s', idx(ei), s.subj_basename, length(s.proc_removeComps), num2str(s.proc_removeComps)) );
        
    end
    
end

if flg == 2, save_desc = 'in_progress'; else, save_desc = 'Default'; end


obj.createResultsCsv(obj.sub, 'postcomps', save_desc);
%obj.createComponentCsv(obj.sub);
mm( sprintf('CSV & mat file: %s\n', obj.htpcfg.csvfile));

%obj.msgout(sprintf('Summary'));
%obj.msgout(sprintf('\nFiles Processed: %d\n', length(objStageStatus)));
%obj.msgout(sprintf('Previously Completed Files: %d\n', prev_files));
%obj.msgout(sprintf('Skipped Files: %d\n', skip_files));
%obj.msgout(sprintf('Other Errors: %d\n', errorchk));

% obj.createResultsCsv(obj.sub, 'postcomps', 'Default');

end

%% disp_welcome
% Display beginning notification for stage 4 of preprocessing to user 
function disp_welcome(mm)
    mm('Starting Stage 4: Post-ICA Component Identification');
end