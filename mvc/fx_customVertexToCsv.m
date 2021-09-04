vertexCsv = bstsub.ImageGridAmp;
timeCsv = bstsub.Time;
vertexCsv = vertexCsv(1:1000, 1:1000);
timeCsv = timeCsv(1:1000)';

csvchan = [];
for chani = 1 : size(vertexCsv,1)
    
    timestrip = vertexCsv(chani, :)';
    chanstrip = repmat(chani, length(timestrip),1);
    
    if chani == 1
        csvchan = [chanstrip timeCsv timestrip];
    else
        csvchan = [csvchan; [chanstrip timeCsv timestrip]];
    end
    
end

%%
parquetwrite('testP2.parquet', array2table(vertexCsv));
writematrix(csvchan, 'test.csv');
csvchan = num2str(csvchan);
fid = fopen('data.txt','w') ;
fprintf(fid,'%s\n',csvchan{:});
fclose(fid) ;
ds = datastore(
M = zeros(1721, 196609);
FileName = 'csvLargeFile.dat';
%To write
vertexCsvT = tall(vertexCsv);
FID = fopen(FileName, 'w');
fwrite(FID, vertexCsv, 'double');
fclose(FID);

tic
fid = fopen('temp4B.txt','wt');
fprintf(fid,csvchan);
fclose(fid);
toc

