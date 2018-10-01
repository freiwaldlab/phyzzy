classdef mlmat < handle
    properties (SetAccess = protected)
        filename
    end
    properties (Access = protected)
        mode
        file_exist
    end
    
    methods
        function obj = mlmat(filename,mode)
            if ~exist('mode','var'), mode = 'r'; end
            if exist('filename','var'), obj.open(filename,mode); end
        end
        function open(obj,filename,mode)
            obj.filename = filename;
            obj.mode = mode;
            obj.file_exist = 2==exist(obj.filename,'file');
        end            
        function close(~), end
        
        function write(obj,val,name)
            if isobject(obj)
                field = fieldnames(val);
                for m=1:length(field)
                    a.(name).(field{m}) = val.(field{m}); %#ok<*STRNU>
                end
            else
                a.(name) = val;
            end
            switch obj.mode
                case 'a'
                    if obj.file_exist
                        save(obj.filename,'-struct','a','-append','-nocompression');
                    else
                        save(obj.filename,'-struct','a','-v7.3','-nocompression');
                        obj.file_exist = true;
                    end
                case 'w'
                    save(obj.filename,'-struct','a','-v7.3','-nocompression');
                    obj.mode = 'a';
                    obj.file_exist = true;
                case 'r'
                    error('This file is read-only!!!');
                otherwise
                    error('Unknown file access mode!!!');
            end
        end
        function val = read(obj,name)
            a = load(obj.filename,name);
            val = a.(name);
        end
        function val = read_trial(obj)
            a = load(obj.filename,'-regexp','^Trial\d+$');
            field = fieldnames(a);
            for m=1:length(field)
                val(m) = a.(sprintf('Trial%d',m)); %#ok<AGROW>
            end
        end
        function val = who(obj)
            val = who('-file',obj.filename);
       end
    end
end
