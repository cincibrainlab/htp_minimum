% bst_clean

function sStat = fx_BstRetrieveGroupStats()
% get computed stats

% identify current study
sStudy = bst_get('Study');

% load names of available stat results
statFiles = {sStudy.Stat.FileName}';

% populate stat result structure
for stat_i = 1 : length( statFiles )
   
   sStatTmp(stat_i) = bst_get('AnyFile', statFiles{stat_i});
   statFileName = fullfile(bst_get('ProtocolInfo').STUDIES, ...
       sStatTmp(stat_i).Stat(stat_i).FileName);
   tmpStatResults = load(statFileName);
   sStat(stat_i) = tmpStatResults;

end

% clean up workspace
clear sStudy stat_i sStatTmp tmpStatResults statFileName statList statFiles

end