function sFiles = process_extract_values_wrapper(sFiles, outfile_csv, subject_number_check)

% extract all values at each vertex
% Process: Extract values: [-1.000s,80.000s] 2-140Hz
sFiles_values = bst_process('CallProcess', ...
    'process_extract_values', sFiles, [], ...
    'timewindow', [-1, 80], ...
    'freqrange',  [2, 90], ...
    'rows',       '', ...
    'isabs',      0, ...
    'avgtime',    1, ...
    'avgrow',     0, ...
    'avgfreq',    0, ...
    'matchrows',  0, ...
    'dim',        2, ...  % Concatenate time (dimension 2)
    'Comment',    '');


% get time frequency values variable
subcsd = in_bst(sFiles_values.FileName); % Freqs: {8×3 cell} TF: [15002×141×8 double] Time: [1×141 double]
sCortex = in_tess_bst(subcsd.SurfaceFile);
subcsd.iAtlas = find(strcmpi({sCortex.Atlas.Name}, 'Desikan-Killiany'));
subcsd.Atlas = sCortex.Atlas(subcsd.iAtlas);
atlas = subcsd.Atlas;

%%
% all subjects
results = {};
%subjectnames = subcsd.History(4:end,3);

 subjectnames = {sFiles(:).SubjectName};

count = 1;
if length(subjectnames) == subject_number_check % hard coded subject number for sanity check
    for i = 1 : length(subjectnames)
        subname = subjectnames{i};
        eegid = extractBetween(sFiles(i).SubjectName, '-N','_p');
        groupid = extractBetween(sFiles(i).SubjectName, 'G','-N');
       
%         % fix name
%             subnametmp = strsplit(subname);
%             subnametmp2 = subnametmp(end);
%             tmpstr = strsplit(subnametmp2{1},'/');
%             tmpstr = tmpstr{1};            
%             bst_file = strsplit(tmpstr,'-');
%             eegid = erase(bst_file{2}, '_postcomp');
%             eegid = regexprep(eegid, 'N','');
%             groupid = bst_file{1};
% 

        for j = 1 : size( atlas.Scouts,2)
            sname =  atlas.Scouts(j).Label;
            sregion =  atlas.Scouts(j).Region;
            svertex =  atlas.Scouts(j).Vertices;
            try
                Freqs = subcsd.Freqs;
                for k = 1 : size(Freqs,1)
                    avgTF = mean(subcsd.TF(svertex,i,k));
                    bandname = Freqs{k,1};
                    results{count,1} = count;
                    results{count,2} = groupid{1};
                    results{count,3} = eegid{1};
                    results{count,4} = sname;
                    results{count,5} = sregion;
                    results{count,6} = bandname;
                    results{count,7} = avgTF;
                    
                    count = count + 1;
                end
                
            catch
                Freqs =  {'chirp', '37,43', 'itpc'};
                avgTF = mean(subcsd.ImageGridAmp(svertex,i));
                bandname = Freqs{1,1};
                results{count,1} = count;
                results{count,2} = groupid{1};
                results{count,3} = eegid{1};
                results{count,4} = sname;
                results{count,5} = sregion;
                results{count,6} = bandname;
                results{count,7} = avgTF;
                
                count = count + 1;
                
            end
           
        end
        
    end
end

% convert (A-m)^2 to (pA-m)^2
results(:,7) = cellfun(@(x) x*1e24, results(:,7), 'un',0);

bstsourcepow = cell2table(results, 'VariableNames', {'id','group',...
    'eegid','label','region','bandname','value'});

writetable(bstsourcepow, fullfile(outfile_csv));

end
