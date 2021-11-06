function [subnames, groupids] = fx_customGetSubNames( sFiles, p, formatType )
% create project specific tags
subnames = {};
groupids = {};

switch formatType
    case "default"
        for i = 1:length(sFiles)
            subnames{i} = sFiles(i).SubjectName;
            groupids{i} = p.getGroupFromSubjectName(subnames{i});
        end
    case "custom1"
        for i = 1:length(sFiles)
            eegid = regexpi(sFiles(i).SubjectName, '(D.+?(?=_postcomp))', 'match');
            groupid = regexpi(sFiles(i).SubjectName, 'Group\d', 'match');
            subnames{i} = eegid{1};
            groupids{i} = groupid{1};
        end
end

end