function s = evp_setfile( s )
% function s = evp_setfile( s,phase_chan )
    % Spectrum
    s.loadDataset('postcomps');

   % s.loadDataset('source');
   % s.createEpochs;
    
%     s.compRemove; % % Remember New datasets does NOT need this line
%     s.setFreqTable([2,3.5;4,7.5;8,10;10.5,12.5;13,30;30.5,55;65,80]); % 7-band SAT
    % s.setFreqTable([.5,3.5; 4,7.5; 8,10; 10.5,12.5; 13,20; 20.5,30; ...
  %                       30.5,58; 62,100; 100.5, 110]); % P1 & CogFlex
%     s.setFreqTable([.5,3.5; 4,7.5; 8,10; 10.5,12.5; 13,20; 20.5,30; ...
%                         30.5,58; 62,110; 130, 170]); % even higher
     s.setFreqTable([2,4; 4,7.5; 8,10; 10.5,12.5; 13,30; 30.5,55;65,80]); % new 11/22
%     s.setFreqTable([1,4; 4,8; 8,13; 13,30; 30,55; 65,100]); % mice
     s.getPntsTable;
     s.generateStoreRoom;
%     s.peakdetect_v1;
     s.bandAverageTrials;
     s.subj_trials = s.EEG.trials;

    opt.gpuOn = 1;
%     s = ram_conn_dwpli( s, opt );
%     if s.EEG.trials >= 15
%         s = modulationIndex_15( s , opt , phase_chan);
%     else
%         str = sprintf(['ERROR: too few data ', s.subj_basename, ...
%             ' number of trials=', num2str(s.EEG.trials), '.']);
%         s.msgout(str,'proc_err');
%     end

%     s.correlation_global_coupling_p1;
%     s.correlation_global_coupling_p2;


% channel
%     for k=1:30
%         chan = ['ch ',num2str(k,'%.2u')];
%         s.complex_wavelet(chan);
%         s.plot_tfmap(chan);
%     end

% region
%     regions = { 16:19, 12:15, [22,25,28], [3,6,9],...
%             [20:21,23:24,26:27,29:30], [1:2,4:5,7:8,10:11], 1:30 };
%     nRegions = length(regions);
%     for k=1:nRegions
%         s.complex_wavelet(regions{k});
%         s.plot_tfmap(k);
%         % only the last region saved
%     end

    s.unloadDataset;
end