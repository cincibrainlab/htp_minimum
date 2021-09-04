function atlas = fx_BstGetDKAtlasFromSurfaceMat()

% retrieve atlas
sCortex = in_tess_bst('@default_subject/tess_cortex_pial_low.mat');
iAtlas = find(strcmpi({sCortex.Atlas.Name}, 'Desikan-Killiany'));
atlas = sCortex.Atlas(iAtlas);

end