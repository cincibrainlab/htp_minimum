function sFiles = fx_BstDeleteFiles( sFiles )

sFiles = bst_process('CallProcess', 'process_delete', sFiles, [], ...
    'target', 1);  % Delete selected files
db_reload_database('current');
end