classdef mvarHandler < handle
    %MVARHANDLER Class to handle group level mvarClass
    
    properties
        g
        resultFileTable;
        subIndex;
        savePath
        savePathAvailable
    end
    
    methods
        function obj = mvarHandler()
            %MVARHANDLER Construct an instance of this class
            %   Pass in GED results
            obj.g = {};
            obj.savePathAvailable = false;
        end
        function [ged,idx] = findSub( obj, eegid )
            for gi = 1 : length( obj.g )
                ged = obj.g{gi};
                subIndex{gi} = ged.s.subj_basename;
            end
            obj.subIndex = subIndex;            
            idx = find(strcmp(eegid, obj.subIndex));
            ged = obj.g{idx};            
        end      
        function obj = setSavePath( obj, pathname )
            obj.savePath = pathname;
            obj.savePathAvailable = true;
        end
        function obj = setPathDb( obj, pathdb )
            ged = obj.g;
            nonEmptyIndex = find(~cellfun(@isempty,ged));
            cellfun(@(x) x.setPathDb(pathdb), ged(nonEmptyIndex), 'UniformOutput', false);
            obj.g = ged;
        end
        function res = getGedArray( obj ), res = obj.g; end
        function obj = loadGedArray( obj, gedArr), obj.g = gedArr; end
        function obj = pushGedObject(obj, gedObj )
            g = obj.getGedArray;
            if isempty( g )
                g{1} = gedObj; 
            else
                g{end+1} = gedObj; 
            end
            obj.g = g;
        end
        function res = getGedObjByIndex(obj, idx)
            res = obj.g{idx};
        end
        function res = getAllGedObj(obj)
            res = obj.g;
        end
        function obj = saveMvarHandler( obj, pathname )
            if obj.savePathAvailable
                for i = 1 : length(obj.g)
                    ged = obj.g{i};
                    ged.clearLastEEG;
                end
                
                formatOut = 'mm_dd_yy_HHMMSS';
                datestr(now,formatOut)
                filename = sprintf('mvarHandlerInstance_%s.mat', ...
                    datestr(now,formatOut));
                fullpath = fullfile(obj.savePath, filename);
                save(fullpath,'obj','-v7.3');
                fprintf('\n\nSUCCESS:\nFile Saved.\n%s\n', fullpath);
            else
                fprintf('\n\nERROR:\nNo save path available.\nSet with setSavePath method.\n\n');
            end
        end
        function obj = loadMvarHandler( obj, filename )
            tmp = load(filename);
            obj = tmp.obj;
        end       
    end
end

