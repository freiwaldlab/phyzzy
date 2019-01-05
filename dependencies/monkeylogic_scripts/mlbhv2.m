classdef mlbhv2 < handle
    properties (SetAccess = protected)
        filename
    end
    properties (Access = protected)
        fid
    end
    properties (Hidden = true)
        var_pos  % [name start end]
        filesize
    end
    
    methods
        function obj = mlbhv2(filename,mode)  % mode: read, write, append
            obj.fid = -1;
            if ~exist('mode','var'), mode = 'r'; end
            if exist('filename','var'), obj.open(filename,mode); end
        end
        function open(obj,filename,mode)
            close(obj);
            if ~exist('mode','var'), mode = 'r'; end
            obj.filename = filename;
            obj.fid = fopen(filename,mode);
            
            obj.var_pos = [];
            fseek(obj.fid,0,1);
            obj.filesize = ftell(obj.fid);
        end            
        function close(obj)
            if -1~=obj.fid, fclose(obj.fid); obj.fid = -1; end
        end
        function delete(obj), close(obj); end
        
        function write(obj,val,name)
            type = class(val);
            if isobject(val), type = 'struct'; end
            switch type
                case 'struct'
                    dim = ndims(val);
                    sz = size(val);
                    field = fieldnames(val);
                    nfield = length(field);
                    fwrite(obj.fid,length(name),'uint64');
                    fwrite(obj.fid,name,'char*1');
                    fwrite(obj.fid,length(type),'uint64');
                    fwrite(obj.fid,type,'char*1');
                    fwrite(obj.fid,dim,'uint64');
                    fwrite(obj.fid,sz,'uint64');
                    fwrite(obj.fid,nfield,'uint64');
                    for m=1:prod(sz)
                        for n=1:nfield, write(obj,val(m).(field{n}),field{n}); end
                    end
                case 'cell'
                    dim = ndims(val);
                    sz = size(val);
                    fwrite(obj.fid,length(name),'uint64');
                    fwrite(obj.fid,name,'char*1');
                    fwrite(obj.fid,length(type),'uint64');
                    fwrite(obj.fid,type,'char*1');
                    fwrite(obj.fid,dim,'uint64');
                    fwrite(obj.fid,sz,'uint64');
                    for m=1:prod(sz), write(obj,val{m},''); end
                case 'function_handle'
                    write_variable(obj,name,func2str(val));
                otherwise
                    write_variable(obj,name,val);
            end
        end
        function val = read(obj,name)
            pos = 0;
            if exist('name','var')
                if ~isempty(obj.var_pos)
                    row = find(strcmp(obj.var_pos(:,1),name));
                    if ~isempty(row)
                        fseek(obj.fid,obj.var_pos{row,2},-1);
                        val = read_variable(obj);
                        return
                    end
                    pos = obj.var_pos{end,3};
                end
            else
                obj.var_pos = [];
            end
            fseek(obj.fid,pos,-1);
            idx = size(obj.var_pos,1);
            while pos < obj.filesize
                try
                    s = ftell(obj.fid);
                    [a,b] = read_variable(obj);
                    idx = idx + 1;
                    obj.var_pos{idx,1} = b;
                    obj.var_pos{idx,2} = s;
                    pos = ftell(obj.fid); obj.var_pos{idx,3} = pos;
                    if exist('name','var')
                        if strcmp(b,name), val = a; return, end
                    else
                        val.(b) = a;
                    end
                catch err
                    warning(err.message);
                    break;
                end
            end
            if ~exist('val','var'), val = []; end
        end
        function val = read_trial(obj)
            if isempty(obj.var_pos)
                pos = 0;
            else
                for m=[obj.var_pos{~cellfun(@isempty,regexp(obj.var_pos(:,1),'^Trial\d+$','once')),2}]
                    fseek(obj.fid,m,-1);
                    [a,b] = read_variable(obj);
                    val(str2double(regexp(b,'\d+','match'))) = a; %#ok<AGROW>
                end
                pos = obj.var_pos{end,3};
            end
            fseek(obj.fid,pos,-1);
            idx = size(obj.var_pos,1);
            while pos < obj.filesize
                try
                    s = ftell(obj.fid);
                    [a,b] = read_variable(obj);
                    idx = idx + 1;
                    obj.var_pos{idx,1} = b;
                    obj.var_pos{idx,2} = s;
                    pos = ftell(obj.fid); obj.var_pos{idx,3} = pos;
                    if ~isempty(regexp(b,'^Trial\d+$','once')), val(str2double(regexp(b,'\d+','match'))) = a; end
                catch err
                    warning(err.message);
                    break;
                end
            end
            if ~exist('val','var'), val = []; end
        end
        function val = who(obj)
            if isempty(obj.var_pos), pos = 0; else pos = obj.var_pos{end,3}; end
            fseek(obj.fid,pos,-1);
            idx = size(obj.var_pos,1);
            while pos < obj.filesize
                try
                    s = ftell(obj.fid);
                    [~,b] = read_variable(obj);
                    idx = idx + 1;
                    obj.var_pos{idx,1} = b;
                    obj.var_pos{idx,2} = s;
                    pos = ftell(obj.fid); obj.var_pos{idx,3} = pos;
                catch err
                    warning(err.message);
                    break;
                end
            end
            if isempty(obj.var_pos), val = []; else val = obj.var_pos(:,1); end
       end
    end
    
    methods (Access = protected)
        function write_variable(obj,name,val)
            dim = ndims(val);
            sz = size(val);
            type = class(val);
            fwrite(obj.fid,length(name),'uint64');
            fwrite(obj.fid,name,'char*1');
            fwrite(obj.fid,length(type),'uint64');
            fwrite(obj.fid,type,'char*1');
            fwrite(obj.fid,dim,'uint64');
            fwrite(obj.fid,sz,'uint64');
            fwrite(obj.fid,val,type);
        end
        function [val,name] = read_variable(obj)
            try
                lname = fread(obj.fid,1,'uint64=>double');
                name = fread(obj.fid,[1 lname],'char*1=>char');
                ltype = fread(obj.fid,1,'uint64=>double');
                type = fread(obj.fid,[1 ltype],'char*1=>char');
                dim = fread(obj.fid,1,'uint64=>double');
                sz = fread(obj.fid,[1 dim],'uint64=>double');
                if strncmp(type,'ml',2), type = 'struct'; end
                switch type
                    case 'struct'
                        nfield = fread(obj.fid,1,'uint64=>double');
                        for m=1:prod(sz)
                            for n=1:nfield, [a,b] = read_variable(obj); val(m).(b) = a; end %#ok<AGROW>
                        end
                        if exist('val','var'), val = reshape(val,sz); else val = []; end
                    case 'cell'
                        val = cell(sz);
                        for m=1:prod(sz), val{m} = read_variable(obj); end
                    otherwise
                        val = reshape(fread(obj.fid,prod(sz),['*' type]),sz);  % fread can handle only a 2-d size arg.
                end
            catch err
                rethrow(err);
            end
        end
    end
end
