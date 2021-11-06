classdef repMakeClass < handle
    %REPMAKE GNU Make for Reproducible Manuscripts
    %   MATLAB Utilities
    
    properties
        Property1
    end
    
    methods
        function obj = repMakeClass()
            %REPMAKE Construct an instance of this class
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
    methods (Static)
        % File Handling: Create target file based on file extension
        function target_file = outFile( prefix, buildPath, typeOfOutput )
            switch typeOfOutput
                case 'MAT'
                    suffix = '.mat';
                case 'CSV'
                    suffix = '.csv';
                case 'PDF'
                    suffix = '.pdf';
                case 'PARQUET'
                    suffix = '.parquet';
                case 'PNG'
                    suffix = '.png';
            end
            target_file = fullfile(buildPath, [prefix suffix]);
        end
    end
    
end

