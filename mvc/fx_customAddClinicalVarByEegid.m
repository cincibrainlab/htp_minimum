function [res idx] = fx_customAddSexVarByEegid( sFiles, csvfile )

[eegids, groupids] = fx_customGetSubNames(sFiles);

clindata = readtable(csvfile);
focustable = [clindata.eegid clindata.group clindata.sex];

for i = 1 : length(eegids)
   current_id = eegids(i);
   matchidx = find(strcmp(current_id, clindata.eegid)); 
   clinvar(i) = clindata.sex(matchidx);
   clinidx{i} = [clindata.group{matchidx} '_'  clindata.sex{matchidx}];
end


for grpi = 1 : length(clinidx)
    
    switch clinidx{grpi}
        case 'FXS_M'
            grpIdxClin(grpi) = 1;
        case 'TDC_M'
            grpIdxClin(grpi) = 2;
        case 'FXS_F'
            grpIdxClin(grpi) = 3;
        case 'TDC_F'
            grpIdxClin(grpi) = 4;
    end
end

res = clinidx;
idx = grpIdxClin;
end