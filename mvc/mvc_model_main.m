% Model Main
% HTP MVC Platform

% index
% id    code      source    description         subjects
% --    ----      ----      ---------------     -----------------
% 01    fxrest    P1        FX resting U54A     71 control, 70 FX
% 02    fxchirp   P1        FX chirp U54A       
% 03    asdrest   apdasd    APD in ASD          20 ASD
% 04    mearest   ucr       WT vs. KO           10 WT, 11 KO
% 05    fxbac     P2        Baclofen/Placebo    17 FX

model_fxrest    = mvcHelperClass('fxrest');
model_fxchirp   = mvcHelperClass('fxchirp');
model_asdrest   = mvcHelperClass('asdrest');
model_mearest   = mvcHelperClass('mearest');
model_fxBac     = mvcHelperClass('fxBac');
model_beats     = mvcHelperClass('beats');

% function dict
% Model
    % assignDirectories(importdir, exportdir)
% Controller
% View
importdir.pc =  
exportdir.pc = 
% Model: FX Resting Dataset (Resting)
model_fxrest    = mvcHelperClass('fxrest');
model_fxrest.assignDirectories(, '/Volumes/extm1/data/fxrest' );
model_fxrest.generateSyspath();
model_fxrest.assignModel('A1911071145_subjTable_P1Stage4.csv');
model_fxrest.createNewFilenameFromBasefile('csv', 'qivars');
model_fxrest.buildModel;



% Model: FX Chirp Dataset (ERP)

% Model: ASD Resting Dataset (Resting)

