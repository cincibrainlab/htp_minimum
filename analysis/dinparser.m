function s = dinparser( s )

% input EEG, din_class, din_target, maxlatency

% output EEG (removed bad trials),
EEG = s.EEG;

%
epoch_n     = length(EEG.epoch); % total no of starting epochs in EEG file
maxlatency = 1050;  % search only DINs up to this time duration
filter = [];    % temp duration filter for each epoch

din_class = {'DIN6', 'DIN7'};   % major DINs to classify trials
din_target = {'DIN9', 'DI11'};  % final DINs to create epochs
din_buttonpress = {'DIN1', 'DIN2', 'DIN3', 'DIN4'};

din_class_n = length(din_class); % array for num of classification DINs
din_target_n = length(din_target); % array for num of target DINs

din_logical_test = 0;

din_valid = zeros( din_class_n, epoch_n );  % array to signify bad trials (1 = marked for removal)
din_bad = [];

din_code = s.eventcfg.nowTrigger;

res = struct(); % results array
res.badtrials = []; % trials removed as not meeting criteria
res.badtrials_n = 0; % number of bad trials removed
% find trials that don't meet criteria
tmpstr = sprintf('[%s]', din_class{:});

str = fprintf('Step 1: %d total trials epoched by %s', epoch_n, s.eventcfg.nowTrigger{1}); s.msgout(str,'step_msg');
for i = 1 : epoch_n
    % examine events before 1050
    filter = (cell2mat(EEG.epoch(i).eventlatency) < maxlatency);
    % check each din_class if present in target window
    % set for two DIN 6 and 7
    
    if i==1, str = sprintf('Step 2: Find valid categories for %s before %d ms.', tmpstr, maxlatency); s.msgout(str, 'step_msg'); end
    
    for j = 1 : din_class_n
        
        epoch_events = EEG.epoch(i).eventtype(filter);
        din_logical_test(j) = any(strcmp( din_class{j},epoch_events ));
        if din_logical_test(j) == 0
            din_valid(j, i) = 0;
        else
            if din_logical_test(j) == 1
                din_valid(j, i) = 1;
            end
        end
    end
    
    % create marked trial for removal that do not contain class
    
    if  any(din_logical_test)
        %  EEG.epoch(i).eventtype(filter)
        din_bad(i) = 0;
    else
        din_bad(i) = 1;
        %[EEG.epoch(i).eventtype; EEG.epoch(i).eventlatency]
        % badtrials{end+1} = [EEG.epoch(i).eventtype;  EEG.epoch(i).eventlatency];
    end
    
    
end

res.badtrials = find(din_bad);
res.badtrials_n = length(find(din_bad));

str = sprintf('\tBad Trial ID: %s\tBad Trial Count: %d' ,mat2str(res.badtrials), res.badtrials_n); s.msgout(str, 'step_complete');
str = sprintf('\tRejecting %d bad trials and creating new dataset.', res.badtrials_n); s.msgout(str, 'step_complete');

% reject bad trials (use EEGLAB function for data integrity)
EEG2 =  pop_rejepoch( EEG, din_bad ,0);
EEG2 = eeg_checkset( EEG2 );

% get index arrays for each din
epoch_n     = length(EEG2.epoch); % total no of starting epochs in EEG file
maxlatency = 1050;  % search only DINs up to this time duration
filter = [];

[a,b,c] = fileparts(EEG.filename); EEG2.setname = [b '_CL' c]; EEG2.filename = EEG2.setname;

str = sprintf('Step 3: Separate remaining %d trials by DIN category (%s) into %s.', epoch_n, tmpstr, EEG2.setname); s.msgout(str, 'step_complete');

for i = 1 : epoch_n     
    filter = (cell2mat(EEG2.epoch(i).eventlatency) < maxlatency);

    % check each din_class if present in target window
    for j = 1 : din_class_n
        %disp(din_class{j});
        epoch_events = EEG2.epoch(i).eventtype(filter);
        

        din_test = [];
        %din_test = any(strcmp( epoch_events, din_class{j}));
        din_old(i,j) = any(strcmp( epoch_events, din_class{j}));
        
        if isequal(i,1) || isequal(i,2)
            if strcmp(din_class{j},'DIN6')
                %din_test = any(strcmp( epoch_events, din_class{j})) && ~any(strcmp(din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)},string([EEG2.epoch(i+1:i+2).eventtype])));
                din_test = any(strcmp( epoch_events, din_class{j})) && ~any(strcmp( epoch_events, din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)})) && ~any(strcmp(din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)},string([EEG2.epoch(i+1:i+2).eventtype])));
            else
                %din_test = any(strcmp( epoch_events, din_class{j})) && ~any(strcmp(din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)},string([EEG2.epoch(i+1).eventtype])));
                din_test = any(strcmp( epoch_events, din_class{j})) && ~any(strcmp( epoch_events, din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)})) && ~any(strcmp(din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)},string([EEG2.epoch(i+1).eventtype])));
            end
        elseif isequal(i,epoch_n) || isequal(i,epoch_n-1)
            if strcmp(din_class{j},'DIN6')
                %din_test = any(strcmp( epoch_events, din_class{j})) && ~any(strcmp(din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)},string([EEG2.epoch(i-2:i-1).eventtype])));
                din_test = any(strcmp( epoch_events, din_class{j})) && ~any(strcmp( epoch_events, din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)})) && ~any(strcmp(din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)},string([EEG2.epoch(i-2:i-1).eventtype])));
            else
               %din_test = any(strcmp( epoch_events, din_class{j})) && ~any(strcmp(din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)},string([EEG2.epoch(i-1).eventtype])));
               din_test = any(strcmp( epoch_events, din_class{j})) && ~any(strcmp( epoch_events, din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)})) && ~any(strcmp(din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)},string([EEG2.epoch(i-1).eventtype])));
            end
        else 
           if strcmp(din_class{j},'DIN6')
%                din_test = any(strcmp( epoch_events, din_class{j})) && ~any(strcmp(din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)},string([EEG2.epoch(i+1:i+2).eventtype]))) && ...
%                    ~any(strcmp(din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)},string([EEG2.epoch(i-2:i-1).eventtype])));
               din_test = any(strcmp( epoch_events, din_class{j})) && ~any(strcmp( epoch_events, din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)})) && ~any(strcmp(din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)},string([EEG2.epoch(i+1:i+2).eventtype]))) && ...
                   ~any(strcmp(din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)},string([EEG2.epoch(i-2:i-1).eventtype])));
           else
%                din_test = any(strcmp( epoch_events, din_class{j})) && ~any(strcmp(din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)},string([EEG2.epoch(i+1).eventtype]))) && ...
%                    ~any(strcmp(din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)},string([EEG2.epoch(i-1).eventtype])));
               din_test = any(strcmp( epoch_events, din_class{j})) && ~any(strcmp( epoch_events, din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)})) && ~any(strcmp(din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)},string([EEG2.epoch(i+1).eventtype]))) && ...
                   ~any(strcmp(din_class{~cellfun(@(x) strcmp(din_class{j},x),din_class)},string([EEG2.epoch(i-1).eventtype]))); 
           end
        end
        din_new(i,j) = din_test; 
        %din_logical_test(j) = any(strcmp( din_class{j},epoch_events ));
        if din_test == 1
            din_valid(j, i) = 1;
        else
            if din_test == 0
                din_valid(j, i) = 0;
            end
        end
        
    end
    
    
end

for i = 1 : size(din_valid,1)
    % create EEG arrays for each DIN
    EEGTEMP = pop_rejepoch( EEG2, ~din_valid(i,:) ,0); % reverse array to select rejections
    EEGDIN(i) = eeg_checkset( EEGTEMP );
    %s.msgout(sprintf('Identified Event %d: [%s] Total Trials: %d', i, din_class{i}, length(EEGDIN(i).epoch)),'step_complete');
    
    
    [a,b,c] = fileparts(EEGDIN(i).filename);
    EEGDIN(i).filename = [a b '_' din_class{i} '_' c];
    EEGDIN(i).setname = EEGDIN(i).filename;
    
end
tmpstr2=sprintf('[%s]', din_target{:});
tmpstr3=sprintf('[%s]', din_buttonpress{:});
str = sprintf('Step 4: Separate each of the %d categories (%s) into success/fail trials based on target (%s)', length(din_class), tmpstr, tmpstr2); s.msgout(str, 'step_complete');

for k = 1 : size(EEGDIN,2)
    % get index arrays for each din
    epoch_n     = length(EEGDIN(k).epoch); % total no of starting epochs in EEG file
    maxlatency_class = 1050;  % search only DINs up to this time duration
    maxlatency_target = 3000;
    filter_class = [];
    filter_target = [];
    din_target_idx = zeros(length(din_target), epoch_n);
    valid_trials = zeros(din_target_n, length(EEGDIN(k).epoch));
    str = sprintf('\tIdentify trials with button (%s) before %d ms', tmpstr3,maxlatency_target); s.msgout(str, 'step_complete');
    for i = 1 : epoch_n
        %disp(i)
        filter_class = (cell2mat(EEGDIN(k).epoch(i).eventlatency) < maxlatency_class); % possible categories
        
        % check each din_class if present in target window
        for j = 1 : din_target_n
            %fprintf('Button Error, Block: %s, Epoch: %d Target: %s\n', din_class{k}, i, din_target{j});
            
            button_presses = [];
            button_no = [];
            
            % events of potential categories, then identify the index of
            % the category
            epoch_events = EEGDIN(k).epoch(i).eventtype(filter_class); epoch_latency = EEGDIN(k).epoch(i).eventlatency(filter_class);
            din_class_idx = find(strcmp( din_class{k}, epoch_events ));
            
            if ~isempty(din_class_idx) % successful identification of category
                din_class_latency = EEGDIN(k).epoch(i).eventlatency{din_class_idx}; % position of the category
                x_target = cell2mat(EEGDIN(k).epoch(i).eventlatency); % temp mat for clarity
                filter_target = ( x_target > din_class_latency & x_target < maxlatency_target); % possible targets
                
                epoch_target = EEGDIN(k).epoch(i).eventtype(filter_target);
                %epoch_events = EEGDIN(k).epoch(i).eventtype(filter_target);
                
                din_target_test = any(strcmp( din_target{j}, epoch_target ));
                target_type = strcmp( din_target{j}, epoch_target );
                
                if din_target_test == 1
                    target_latency_temp = cell2mat(EEGDIN(k).epoch(i).eventlatency(filter_target));
                    target_latency = target_latency_temp(target_type);
                    filter_target2 = ( x_target > din_class_latency & x_target < target_latency); % possible targets
                    
                    din_target_tmp_type = EEGDIN(k).epoch(i).eventtype(filter_target2);
                    din_target_tmp_latency = EEGDIN(k).epoch(i).eventlatency(filter_target);
                    button_presses = ismember( din_target_tmp_type, din_buttonpress);
                    button_no = length(any( button_presses ));
                    
                    if button_no > 1
                        filter_disp = (cell2mat(EEGDIN(k).epoch(i).eventlatency) < maxlatency_target); % possible categories
                        
                        str = sprintf('[%s]', EEGDIN(k).epoch(i).eventtype{filter_disp});
                        str2 = sprintf('Trial %d: %s', i, str);
                        str_button = sprintf('[%s]', din_target_tmp_type{button_presses});
                        str3 =  sprintf('Error: Too many button presses. Block: %s, Epoch: %d Target: %s Button: %s\n', din_class{k}, i, din_target{j}, str_button);
                     %   s.msgout(str2,'proc_complete');
                    %    s.msgout(str3,'proc_complete');
                    end
                    if button_no == 0
                        s.msgout('No button presses.','proc_complete');
                        filter_disp = (cell2mat(EEGDIN(k).epoch(i).eventlatency) < maxlatency_target); % possible categories
                        
                        str = sprintf('[%s]', EEGDIN(k).epoch(i).eventtype{filter_disp});
                        str2 = sprintf('Trial %d: %s', i, str);
                        %str_button = sprintf('[%s]', din_target_tmp_type{button_presses});
                        str3 =  sprintf('Error: NO. Block: %s, Epoch: %d Target: %s\n', din_class{k}, i, din_target{j});
                     %   s.msgout(str2,'proc_complete');
                     %   s.msgout(str3,'proc_complete');
                    end
                    if button_no == 1
                        filter_disp = (cell2mat(EEGDIN(k).epoch(i).eventlatency) < maxlatency_target); % possible categories
                        
                        valid_trials(j,i) = 1;
                        str = sprintf('[%s]', EEGDIN(k).epoch(i).eventtype{filter_disp});
                        str2 = sprintf('Trial %d: %s', i, str);
                        str_button = sprintf('[%s]', din_target_tmp_type{button_presses});
                        str3 =  sprintf('Valid Trial. Block: %s, Epoch: %d Target: %s Button: %s', din_class{k}, i, din_target{j}, str_button);
                     %   s.msgout(str2,'proc_complete');
                     %   s.msgout(str3,'proc_complete');
                     %   s.msgout('One button pressed.','proc_complete');
                    end
                    
                    
                    
                else
                    filter_disp = (cell2mat(EEGDIN(k).epoch(i).eventlatency) < maxlatency_target); % possible categories

                    str = sprintf('[%s]', EEGDIN(k).epoch(i).eventtype{filter_disp});
                    str2 = sprintf('Trial %d: %s', i, str);
                  %  s.msgout(str2,'proc_complete');
                    str4 =  sprintf('No valid target. Block: %s, Epoch: %d, Expected Target: %s', din_class{k}, i, din_target{j} );
                   % s.msgout(str4,'proc_complete');
                    
                end
            end
            
            
            
            
        end
        
        
    end
    final_trials{k,1} = EEGDIN(k).setname;
    final_trials{k,2} = valid_trials;
    
end

for i = 1 : size(final_trials,1)
    
    EEG_setname_temp = EEGDIN(i).setname;
    
 
    s.eventcfg.din_code = '';
     
    for j = 1 : size(final_trials{i,2}, 1)
    
        din_code = strcat(s.eventcfg.nowTrigger, '_', din_class{i}, '_',  din_target{j});
        s.eventcfg.din_code = din_code;
        s.eventcfg.details(end+1,1) = din_code;
        
        tmp_fn = sprintf('%s', din_target{j});
        
        [a,b,c] = fileparts(EEG_setname_temp);
        EEGDIN(i).filename = [a b tmp_fn '_' c];
        EEGDIN(i).setname = EEGDIN(i).filename;
        
        tmp_rej = final_trials{i,2};
        str = sprintf('DIN: %s, Total: %d', din_target{j}, length(find(tmp_rej(j,:))));
        s.msgout(str,'proc_complete');
        EEGTMP = pop_rejepoch( EEGDIN(i), ~tmp_rej(j,:) ,0);
        EEGTMP = eeg_checkset( EEGTMP );
        
        
        s.EEG = EEGTMP;
        s = s.compRemove; % remove components
        s.EEG.etc.din_code = din_code;
        s.eventcfg.nowTrigger = din_target(j);
        s.eventcfg.stage = 'C';
        s.createEpochsERP;
        

        %s.storeDataset( s.EEG, s.pathdb.postcomps, s.subj_subfolder, s.EEG.filename);
        
    end

    EEGDIN(i).setname = EEG_setname_temp;

end


end






















