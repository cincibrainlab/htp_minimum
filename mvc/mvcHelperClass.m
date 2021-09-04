classdef mvcHelperClass < handle & htpPortableClass
    %htp mvc Framework Helper
    
    properties
        % model
        projectcode   % dataset identifier
        importdir     % htp data directory
        exportdir     % export directory for analytics
        syspath       % completed path structure
        modelcsv      % csv model from analysis directory
        modelmat      % associated object file
        basefiles     % basefile templates
        outfile       % output cellvector of generated files
        
    end
    
    methods
        function obj = mvcHelperClass(projectcode)
            %mvcHelperClass Construct an instance of this class
            obj.projectcode = projectcode;
            obj.syspath = struct();
            obj.outfile = {'Type', 'Description', 'DateTime','Location'};
        end
        
        function obj = assignDirectories(obj, importdir, exportdir)
            %mvcHelperClass Construct an instance of this class
            obj.importdir   = importdir;
            obj.exportdir   = exportdir;           
        end
    end
    
    % Model
    methods
        
        function syspath = ...
                generateSyspath(obj)
            
            makePath = @(basedir, subfolder) fullfile(basedir, subfolder);
            
            syspath.base        = obj.importdir;
            syspath.base_cloud  = obj.exportdir;
            syspath.dataset     = makePath(syspath.base, 'A00_ANALYSIS'); % cleaned data
            syspath.perms       = makePath(syspath.base, 'perm');
            syspath.figures     = makePath(syspath.base, 'figs');
            syspath.mat         = makePath(syspath.base, 'mat');
            syspath.R           = makePath(syspath.base_cloud, 'R/rawdata');
            syspath.script      = makePath(syspath.base_cloud, 'MATLAB');
            syspath.class       = makePath(syspath.base_cloud, 'MATLAB/classes');
            syspath.external    = makePath(syspath.base_cloud, 'MATLAB/external');
            
            obj.syspath = syspath;
            obj.generateBasefileTemplates();
            
        end
        
        function obj = assignModel( obj, csvfile )
        % model defined from CSV file in A00_Analysis folder
        
        obj.modelcsv = fullfile(obj.syspath.dataset, csvfile);
        obj.modelmat = obj.csvfilename2matfilename(obj.modelcsv);
        
        end
        
        function obj = generateBasefileTemplates( obj )
            obj.basefiles = struct();
            obj.basefiles.csv = fullfile(obj.syspath.R, [obj.projectcode '_csv_TBD_TIME.csv']);
            obj.basefiles.fig = fullfile(obj.syspath.R, [obj.projectcode '_figfile_TBD_TIME.fig']);
            obj.basefiles.svg = fullfile(obj.syspath.R, [obj.projectcode '_figfile_TBD_TIME.svg']);
            obj.basefiles.pdf = fullfile(obj.syspath.R, [obj.projectcode '_figfile_TBD_TIME.pdf']);
            obj.basefiles.m = fullfile(obj.syspath.R, [obj.projectcode '_mfile_TBD_TIME.m']);
            obj.basefiles.mat = fullfile(obj.syspath.R, [obj.projectcode '_matfile_TBD_TIME.mat']);
            
        end    
        
        function filename = createNewFilenameFromBasefile( obj, type, description)
            datecode = datestr(datetime('now'), 'mmddyyhhMM');
            filename = strrep(obj.basefiles.(type), 'TBD', description); 
            filename = strrep(filename, 'TIME', datecode); 
            obj.outfile{end+1,1} = type;
            obj.outfile{end, 2} = description;
            obj.outfile{end, 3} = datecode;
            obj.outfile{end, 4} = filename; 

        end
        
        function obj = generateFilepath(obj)
            obj.outfile.csv
            
            outfile.csv.prep    = makePath(syspath.R, '00_gedB_QI_preprocessing.csv');
            outfile.csv.trials  = makePath(syspath.R, '01_gedB_QI_validtrials.csv');
            
            
        end
        
        function obj = buildModel(obj)
            obj.importDataObjects(obj.modelcsv, obj.modelmat, obj.syspath.dataset);
            obj.updateBasePaths( obj.syspath.dataset );
        end
    end
    
    % Controller
    methods
    end
    
    % View
    methods
    end
    
    % Static
    methods (Static)
        
        function res = makeGenericPath(basedir, subfolder)
            res = fullfile(basedir, subfolder);
        end
        
        function res = makeCodeString( code, str)
            res = sprintf("%s_%s", code, str);
        end
        
        function matfile = csvfilename2matfilename( csvfilename )
            [path, name, ~] = fileparts(csvfilename);
            matfile = fullfile(path, [name '.mat']);
        end
    end
end

