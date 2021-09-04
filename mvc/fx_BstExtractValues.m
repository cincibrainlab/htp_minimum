function resulttable = fx_BstExtractValues( sFiles )

[subnames, groupids] = fx_customGetSubNames(  sFiles );

% get TF values
sValues = bst_process('CallProcess', ...
    'process_extract_values', sFiles, []);
sTf = in_bst(sValues.FileName);

% retrieve atlas
sCortex = in_tess_bst(sTf.SurfaceFile);
sTf.iAtlas = find(strcmpi({sCortex.Atlas.Name}, 'Desikan-Killiany'));
sTf.Atlas = sCortex.Atlas(sTf.iAtlas);
atlas = sTf.Atlas;

% all subjects
cnt = 1;
results = cell(length(subnames)*size( atlas.Scouts,2)*size(sTf.Freqs,1),7);
for i = 1 : length(subnames)
    eegid = subnames{i};
    groupid = groupids{i};
    
    for j = 1 : size( atlas.Scouts,2)
        sname =  atlas.Scouts(j).Label;
        sregion =  atlas.Scouts(j).Region;
        svertex =  atlas.Scouts(j).Vertices;
        Freqs = sTf.Freqs;
        for k = 1 : size(Freqs,1)
            avgTF = mean(sTf.TF(svertex,i,k));
            bandname = Freqs{k,1};
            results{cnt,1} = cnt;
            results{cnt,2} = groupid;
            results{cnt,3} = eegid;
            results{cnt,4} = sname;
            results{cnt,5} = sregion;
            results{cnt,6} = bandname;
            results{cnt,7} = avgTF;
            cnt = cnt + 1;
        end
    end
    
end
resulttable = results;

end