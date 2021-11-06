%% Stage Direct
%%% Usage
% obj = preprocess_direct( obj, cfg )
%
%%% Parameters
%
% * INPUTS: obj, cfg
%  
% * OUTPUTS: obj
%
% The optional input, obj if the function is not self-invoked, is the 
% htpPreprocessMaster object.  
% The output, obj if the function is not self-invoked, is the updated 
% htpPreprocessMaster object with information updated regarding the status 
% and details of stage 1 completion of preprocessing.
%
%%% Copyright and Contact Information
% Copyright � 2020 Cincinnati Children's (Pedapati Lab)
%
% This file is part of High Throughput Pipeline (HTP)
%
% See https://bitbucket.org/eped1745/htp_stable/src/master/
%
% Contact: ernest.pedapati@cchmc.org

%%  preprocess_direct
% entry point function to import completed preprocessed files
% into a htp object for further processing. The use case for this
% function was importing completed auditory ERP trial level data
% to create a htpPreprocessing object for source localization.
% File was based off of preprocess_stage1.m. "direct" stands for
% direct import.
% Configure logs to produce needed notifications, warning, and errors to
% user during stage direct preprocessing.
% 
% Setup stage to start by configuring datasets (channel, user, title, and resample hz),
% loading all RAW files, performing necessary filtering and resampling, 
% and creating a corresponding result csv for each file. 

function obj = preprocess_direct( obj, cfg )

% logging functions
obj.reset_msg_log;                      
[mc, mm, mw] = obj.tools_log;           
% mc(getmsg(obj, 7));       

% add direct import parameters to htpcfg
obj.htpcfg.direct.cfg = cfg;

% stage assigment for CSV
stage_last = 'raw';
stage_next = 'postcomps';

% manual configuration from helper script
% data directory per HTP specs
obj.htpcfg.basePath = cfg.basePath;

% key parameters for direct file import
% channel structure to apply to file
% add import configuration to cfg.htpEegSystems.xml
obj.htpcfg.chanNow = obj.htpcfg.chaninfo(1);
obj.htpcfg.chanNow.net_filter = '*.set';
obj.setOptNow(repmat(1,[1 25]));
% sample rate
obj.htpcfg.optnow.Stage1_Resample =  cfg.srate;
obj.htpcfg.chanNow.net_name = 'SET'
% import description details
obj.htpcfg.user = 'DIRECT';
obj.study_title = cfg.study_title;

% get filelist and subfolders
obj.import_filelist;
fnlist     = obj.fnlist;                
flength    = length(fnlist);
obj.htpcfg.xmax_arr = [];

obj.assignConfig2Subjects;

arrayfun( @(s) setAssets( obj, s ), obj.sub, 'uni', 0);

totalsubs = length( obj.sub );
sub = obj.sub;


for i = 1 : totalsubs
    s = sub(i);
    importData( obj, s );
end

mc(getmsg(obj, 6));                     

obj.sub = sub;
clear sub;

try
    obj.createResultsCsv(obj.sub, stage_next,'Default');
catch
    htps1_backup = createbackup( obj, stage_next );
    assignin('base', 'htps1_backup', htps1_backup);
    try
        save(fullfile(obj.htpcfg.pathdb.analysis, ['backup_' obj.htpcfg.timetag2 '.mat']), 'htps1_backup');
    catch
        obj.msgout('Check folder or reload to create directories.', 'step_error');
    end
end

end

%%  createbackup
% Creates backup of object to be utilized for reprocessing if the initial
% processing does not successfully complete
function backup = createbackup( obj, stage_next )

backup.obj = obj;
backup.stage_next = stage_next;

end

%%  getmsg
% Produces the stage related information that may be of interest to the user
% The command console will display this appropriate information.
function str = getmsg(obj, index )
sf = @sprintf;

opt = obj.formatOptions;
try
    msg{1} = sf('version: htp2020v1\nStarting Stage 1: Import Raw Files...');
    msg{2} = sf('Assign Electrode System: %s', obj.htpcfg.chanNow.net_name);
    msg{3} = sf('User: %s', obj.htpcfg.user);
    msg{4} = sf('Study Title: %s', obj.study_title);
    msg{5} = sf('Set Resample Rate: %d', opt.srate);
    msg{6} = 'Convert RAW data to EEGLAB datasets.';
catch
    msg{7} = sf('Assuming commandline mode, skipping GUI options.');
end

str = msg{index};

end

%%  setAssets
% Configures the information needed for stage 1 preprocessing epr subject
% such as channel files, user, study title, resampling rate, and gui-set
% options
function s = setAssets(obj, s)
[mc, mm, mw] = obj.tools_log;           

opt = obj.formatOptions;
if isempty(fieldnames(opt))
    mw(sprintf('WARNING: Assuming direct (commandline) mode, skipping GUI options.'));
    opt.srate = obj.htpcfg.optnow.Stage1_Resample;
end
    
s.setElectrodeSystem( obj.htpcfg.chanNow );     
s.setUser(obj.htpcfg.user);                     
s.changeStudyTitle(obj.study_title);            
s.setResampleRate( opt.srate );                 
s.setopt( opt );                               

end

%%  importData
% imports the subject data and performs the necessary steps of filtering,
% resampling, and downgrading, as needed.  Then stores updated data and
% configurations for later use in preprocessing.
function s = importData( obj, s )
stage_next = 'postcomps';
opt = obj.formatOptions;

s.getRawData;                                    % import data (biosig)

s.storeDataset( s.EEG, ...
            obj.htpcfg.pathdb.(stage_next), ...
            s.subj_subfolder, ...
            s.filename.( stage_next ) );
        
s.unloadDataset;

s.outputRow( stage_next );
obj.update_htpcfg( s );

end