classdef csvOutClass < handle
    %csvOutClass Helpful object for generating output files
      
    properties
        pathdb
        filename
        filedb
        datatbl
        columnNames
        numRows
        
        
    end
    
    methods
        function obj = csvOutClass()            
            obj.pathdb.R = '';
            obj.pathdb.script = '';
            obj.filedb = struct();
        end
        function res = getNumRows( obj )
            res = obj.numRows;
        end
        function obj = setNumRows(obj, num)
            obj.numRows = num;
        end
        function res = createCellArr(obj)
            if isempty(obj.numRows) || isempty(obj.columnNames)
                disp('Need Column Variable names and/or number of rows.');
            else
            res = cell(obj.numRows, size(obj.columnNames,2));
            % obj.datatbl = table(size(res));
            end
            
        end
        function obj = storeTableFromCell(obj, cellarr)
            obj.datatbl = cell2table(cellarr, 'VariableNames', obj.columnNames);
        end
        function res = getDataTable(obj)
            res = obj.datatbl;
        end
        function obj = setFileName(obj, fileName)
            obj.filename=fullfile(fileName);
        end
        function obj = setDir(obj, pathType, pathString)
            if nargin < 2
                disp('Send fieldname as character vector.');
                disp(fieldnames(obj.pathdb));
            else
                obj.pathdb.(pathType) = fullfile(pathString);
            end
            
        end
        function obj = setDataTbl(obj, tab )
            obj.datatbl = tab;
        end  
        function obj = setColumnNames(obj, arr)
            obj.columnNames = arr;
        end
          function res = getColumnNames(obj)
            res = obj.columnNames;
          end
          function obj = exportDataTable( obj, pathType )
              if nargin < 2
                  disp('Send fieldname as character vector.');
                  disp(fieldnames(obj.pathdb));
              else
                  obj.filedb.(pathType) = fullfile(obj.pathdb.(pathType), ...
                      obj.filename);
                  writetable(obj.datatbl, obj.filedb.(pathType));
                  disp(['Writing ' obj.filedb.(pathType)]); 
              end
          end
              
    end
end

