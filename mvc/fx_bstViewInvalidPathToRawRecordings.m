function currentpath = fx_bstViewInvalidPathToRawRecordings(sFilesRecordings)
% get current filepath to raw data
[subnames, groupids] = fx_customGetSubNames( sFilesRecordings );

statFileName = {};

for subi = 1 : length(sFilesRecordings)
    
    sMatrix = in_bst_data( sFilesRecordings(subi).FileName);
    
    sMatrix.F.filename(strfind(sMatrix.F.filename,'\'))='/';
    
    [invalidpath, matfile, ext] = fileparts(sMatrix.F.filename);
    
    currentpath{subi} =  sMatrix.F.filename;
    
end

currentpath = currentpath';
end