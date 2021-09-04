%% Chirp High Resolution Cortical Investigation

sStudy = bst_get('Study');
sProtocol = bst_get('ProtocolInfo');
sSubjects =  bst_get('ProtocolSubjects');
sStudyList = bst_get('ProtocolStudies'); 
atlas = fx_BstGetDKAtlasFromSurfaceMat;
sCortex = in_tess_bst('@default_subject/tess_cortex_pial_low.mat');

% mean source file to modify
'Group_analysis/sensorychirp/results_MN_EEG_210620_1505_abs_avg.mat'


% get pointers to raw recordings
sFilesRecordings = bst_process('CallProcess', ...
    'process_select_files_data', [], []);

% create epochs
% Process: Import MEG/EEG: Time

sFilesRecordings = fx_bstSplitRecording( sFilesRecordings, 3.252, 'sensorychirp', [] );

% select trial sources
sFilesSource = fx_getBstSourceFiles();

% 1. Relative Power
sFilesSource2 = fx_getBstSelectByTag(sFilesSource, ...
    'sensorychirp', 'select', 10696);




tVec = ((1:size(trialData,3))/500 ) - .502;  
extractVector =  find(tVec >.5 & tVec < 1.250);


[subnames groupids] = fx_customGetSubNames(sFilesSource);

for i = 1 : length(sFilesSource)

tmpData = in_bst(sFilesSource(i).FileName);
trialData = reshape(tmpData.ImageGridAmp, 15002, [], 1626);
selectData = trialData(:,:,extractVector);

resultCell{i} = selectData;

end

allErsp1 = cell(length(sFilesSource), size(tmpData.ImageGridAmp,1)); 
allitc = cell(length(sFilesSource), size(tmpData.ImageGridAmp,1));

for i = 1 : 1 %length(sFilesSource)
    
    tmpData = in_bst(sFilesSource(i).FileName);
    for j = 1 : size(tmpData.ImageGridAmp,1)
        j
        trials = 176;
        [ersp1,itc,n2,t_s,f_s,erspboot,itcboot,tfdata]= ...
            newtimef(tmpData.ImageGridAmp(sensor,:), ...
            1626,[750 1000],500,[1 10],...
            'winsize',100,'nfreqs',10,...
            'freqs',[35 45],...
            'plotersp','off','plotitc','off',...
            'verbose','off','baseline',NaN,'timesout',250);
        
        allErsp1{i,j} = ersp1;
        allitc{i,j} = itc;
    end
    
    
    
end

for i=1:400;
    rcrits(i,1)=sqrt(-(1/i)*log(.5));
end


ITC1=(abs(itc))-rcrits(trials);
ITC1=(abs(itc))-rcrits(trials);

figure; imagesc(t_s,f_s,ersp1), axis xy;
 figure; imagesc(t_s,f_s,ITC1), axis xy;


% sanity check
alldata = permute(trialData, [1 3 2]);
plot(tVec', squeeze(trialData(1,1,:)));

% convert to matlab code for speed
ts = tmpData.ImageGridAmp(sensor,:);
ts_trials = reshape(ts, [], 176);
size(ts_trials)
tVec = ((1:size(trialData,3))/500 ) - .502;  

plot(tVec, ts())

286176/1626

