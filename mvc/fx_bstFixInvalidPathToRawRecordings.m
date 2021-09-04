function fx_bstFixInvalidPathToRawRecordings(sFilesRecordings, correctPath)
% get current filepath to raw data
[subnames, groupids] = fx_customGetSubNames( sFilesRecordings );

statFileName = {};

for subi = 1 : length(sFilesRecordings)
    
    sMatrix = in_bst_data( sFilesRecordings(subi).FileName);
    
    sMatrix.F.filename(strfind(sMatrix.F.filename,'\'))='/';
    
    [invalidPath, matfile, ext] = fileparts(sMatrix.F.filename);
    
    statFileNameNew{subi} = fullfile(correctPath, groupids{subi}, [matfile ext]);
    
    sMatrix.F.filename = statFileNameNew{subi};
    bst_save(fullfile(bst_get('ProtocolInfo').STUDIES, ...
        sFilesRecordings(subi).FileName), sMatrix, 'v7');
    
end
end