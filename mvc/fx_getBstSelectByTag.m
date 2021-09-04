function sFilesSel = fx_getBstSelectByTag( sFiles, tag, selectType, expectedNumber )
% 'search' = {'Search the file paths', 'Search the file names', 'Search the names of the parent file'};
% 'select' = {'Select only the files with the tag', 'Ignore the files with the tag'};

switch selectType
    case 'select'
        select = 1;
    case 'deselect'
        select = 2;
    otherwise
        select = 1;
end

sFilesSel = bst_process('CallProcess', 'process_select_tag', ...
    sFiles, [], ...
    'tag',   tag, ...
    'search', 2, ...
    'select', select);

str = sprintf('Number of files selected: %d', numel(sFilesSel));
disp(str);
disp(unique({sFilesSel.Comment}'));

assert(expectedNumber == length(sFilesSel),'Expected Subject Number does not match Result');
end
