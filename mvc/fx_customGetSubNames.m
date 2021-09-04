function [subnames groupids eegid_base] = fx_customGetSubNames( sFiles )
% create project specific tags
subnames = {};
groupids = {};

for i = 1:length(sFiles)
    eegid = extractBetween(sFiles(i).SubjectName, '-N','_p');
    groupid = extractBetween(sFiles(i).SubjectName, 'G','-N');
    eegid_basetmp = extractBetween(eegid, 'D','_c', 'Boundaries', 'Inclusive');
    eegid_basetmp = extractBetween(eegid_basetmp, 'D','_', 'Boundaries', 'Inclusive');
    eegid_basetmp{1} = [eegid_basetmp{1} 'rest'];
    subnames{i} = eegid{1};
    groupids{i} = groupid{1};
    eegid_base{i} = eegid_basetmp{1};

end

end