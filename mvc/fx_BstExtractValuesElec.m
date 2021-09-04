function [resulttable VariableNames] = fx_BstExtractValuesElec( sFiles )

% multi-spectral power export
[subnames, groupids] = fx_customGetSubNames(  sFiles );

% get TF values
sValues = bst_process('CallProcess', ...
    'process_extract_values', sFiles, []);
sTf = in_bst(sValues.FileName);

% all subjects
cnt = 1;
results = cell(length(subnames)*length( sTf.RowNames)*size(sTf.Freqs,1),7);

% convert units for absolute power
if contains(sFiles(1).Comment, 'rel')
    multiplier = 1;
else
    multiplier = 1e12;
end

if contains(sFiles(1).Comment, 'FFT')
    Freqs = sTf.Freqs;
    [ d, ix ] = min( abs(  sTf.Freqs-2 ) );
    [ d2, ix2 ] = min( abs(  sTf.Freqs-80 ) );
    totalRows = {};
    for i = 1 : length(subnames)
        eegid = subnames(i);
        groupid = groupids(i);
 
        subTFrow = num2cell(squeeze(sTf.TF(:,i,ix:ix2)).* multiplier); % picoAmps
        subTFrow = [repmat(eegid,size(subTFrow,1),1) ...
            repmat(groupid, size(subTFrow,1),1) ...
            sTf.RowNames ...
            subTFrow];
        if i == 1
            totalRows = subTFrow;
        else
            totalRows = vertcat(totalRows,subTFrow);
        end
    end
    resulttable = totalRows;
    VariableNames = [{'eegid'}, {'groupid'}, {'electrode'}, ...
            cellfun(@(x) sprintf('F%2.1f', x),  num2cell(Freqs(ix:ix2)), 'uni', 0)];
else
    for i = 1 : length(subnames)
        eegid = subnames{i};
        groupid = groupids{i};
        
        for j = 1 : length( sTf.RowNames)
            electrodeid =  sTf.RowNames{j};
            
            Freqs = sTf.Freqs;
            for k = 1 : length(Freqs)
                avgTF = mean(sTf.TF(j,i,k));
                bandname = Freqs{k,1};
                results{cnt,1} = cnt;
                results{cnt,2} = groupid;
                results{cnt,3} = eegid;
                results{cnt,4} = sFiles(i).Comment;
                results{cnt,5} = electrodeid;
                results{cnt,6} = bandname;
                results{cnt,7} = avgTF;
                cnt = cnt + 1;
            end
        end
        
    end
    resulttable = results;
    VariableNames = {'id','group',...
    'eegid','method','electrode', 'bandname','value'};
    
end
end