classdef mlhdf5 < handle
    properties (SetAccess = protected)
        filename
    end
    properties (Access = protected)
        fid
        opened
        gcpl
        
        int8
        uint8
        int16
        uint16
        int32
        uint32
        single
        double
    end
    properties (Constant)
        root = '/ML/'
    end
    
    methods
        function obj = mlhdf5(filename,mode)  % mode: r, w, a
            obj.opened = false;
            if ~exist('mode','var'), mode = 'r'; end
            if exist('filename','var'), obj.open(filename,mode); end
        end
        function open(obj,filename,mode)
            close(obj);
            if ~exist('mode','var'), mode = 'r'; end
            obj.filename = filename;

            tracked = H5ML.get_constant_value('H5P_CRT_ORDER_TRACKED');
            indexed = H5ML.get_constant_value('H5P_CRT_ORDER_INDEXED');
            order = bitor(tracked,indexed);
            obj.gcpl = H5P.create('H5P_GROUP_CREATE');
            H5P.set_link_creation_order(obj.gcpl,order);

            switch lower(mode)
                case {'r','read'}
                    if 2~=exist(filename,'file'), error('File not found'); end
                    obj.fid = H5F.open(filename,'H5F_ACC_RDONLY','H5P_DEFAULT');
                case {'w','write'}
                    obj.fid = H5F.create(filename,'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');
                    gid = H5G.create(obj.fid,obj.root,'H5P_DEFAULT',obj.gcpl,'H5P_DEFAULT');
                    H5G.close(gid);
                case {'a','append'}
                    if 2==exist(filename,'file')
                        obj.fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
                    else
                        obj.fid = H5F.create(filename,'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');
                        gid = H5G.create(obj.fid,obj.root,'H5P_DEFAULT',obj.gcpl,'H5P_DEFAULT');
                        H5G.close(gid);
                    end
            end
            obj.opened = true;

            obj.int8 = H5T.copy('H5T_NATIVE_CHAR');
            obj.uint8 = H5T.copy('H5T_NATIVE_UCHAR');
            obj.int16 = H5T.copy('H5T_NATIVE_SHORT');
            obj.uint16 = H5T.copy('H5T_NATIVE_USHORT');
            obj.int32 = H5T.copy('H5T_NATIVE_INT');
            obj.uint32 = H5T.copy('H5T_NATIVE_UINT');
            obj.single = H5T.copy('H5T_NATIVE_FLOAT');
            obj.double = H5T.copy('H5T_NATIVE_DOUBLE');
        end            
        function close(obj)
            if ~obj.opened, return, end
            H5T.close(obj.int8);
            H5T.close(obj.uint8);
            H5T.close(obj.int16);
            H5T.close(obj.uint16);
            H5T.close(obj.int32);
            H5T.close(obj.uint32);
            H5T.close(obj.single);
            H5T.close(obj.double);

            H5P.close(obj.gcpl);
            H5F.close(obj.fid);
            obj.opened = false;
        end
        function delete(obj), close(obj); end
        
        function write(obj,val,name)
            location = [obj.root name];
            type = class(val);
            if isobject(val), type = 'struct'; end
            sz = size(val);
            count = prod(sz);
            switch type
                case 'char', write_char(obj,val,location);
                case 'logical', write_logical(obj,val,location);
                case {'int8','uint8','int16','uint16','int32','uint32','single','double'}, write_numeric(obj,val,location);
                case 'struct'
                    field = fieldnames(val);
                    gid = H5G.create(obj.fid,location,'H5P_DEFAULT',obj.gcpl,'H5P_DEFAULT');
                    write_char_attribute(obj,gid,'class',type);
                    write_numeric_attribute(obj,gid,'size',size(val));
                    write_numeric_attribute(obj,gid,'nfield',length(field));
                    H5G.close(gid);
                    switch count
                        case 1, for m=1:length(field), if ~isempty(val), write(obj,val.(field{m}),[name '/' field{m}]); end, end
                        otherwise, for m=1:count, write(obj,val(m),[name '/' num2str(m)]); end
                    end
                case 'cell'
                    gid = H5G.create(obj.fid,location,'H5P_DEFAULT',obj.gcpl,'H5P_DEFAULT');
                    write_char_attribute(obj,gid,'class',type);
                    write_numeric_attribute(obj,gid,'size',sz);
                    H5G.close(gid);
                    for m=1:prod(sz), write(obj,val{m},[name '/' num2str(m)]); end
                case 'function_handle', write_char(obj,func2str(val),location);
                otherwise
                    error('The variable class, ''%s'', is unknown',type);
            end
        end
        function val = read(obj,name)
            if exist('name','var')
                try val = read_variable(obj,[obj.root name]); catch, val = []; end
            else
                names = who(obj);
                for m=1:length(names)
                    val.(names{m}) = read_variable(obj,[obj.root names{m}]);
                end
            end
        end
        function val = read_trial(obj)
            names = who(obj);
            idx = 0;
            for m=1:length(names)
                if ~isempty(regexp(names{m},'Trial\d+','once'))
                	idx = idx + 1;
                    val(idx) = read_variable(obj,[obj.root names{m}]); %#ok<AGROW>
                end
            end
        end
        function val = who(obj,location)
            if ~exist('location','var'), location = obj.root; end
            group_id = H5G.open(obj.fid,location);
            val = {};
            H5L.iterate(group_id,'H5_INDEX_CRT_ORDER','H5_ITER_INC',0,@iterfunc,0);
            
            function [status,opdata_out] = iterfunc(~,name,~)
                status = []; opdata_out = []; val{end+1} = name;
            end
        end
    end
    
    methods (Access = protected)
        function write_char(obj,val,location)
            type_id = H5T.copy('H5T_C_S1');
            H5T.set_size(type_id,'H5T_VARIABLE');
            space_id = H5S.create('H5S_SCALAR');
            dset_id = H5D.create(obj.fid,location,type_id,space_id,'H5P_DEFAULT');
            H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',{val});
            write_char_attribute(obj,dset_id,'class','char');
            write_numeric_attribute(obj,dset_id,'size',size(val));
            H5D.close(dset_id);
            H5S.close(space_id);
            H5T.close(type_id);
        end
        function write_logical(obj,val,location)
            sz = size(val);
            space_id = H5S.create_simple(ndims(val),fliplr(sz),[]);
            dset_id = H5D.create(obj.fid,location,obj.uint8,space_id,'H5P_DEFAULT');
            H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',uint8(val)); %#ok<CPROPLC>
            write_char_attribute(obj,dset_id,'class','logical');
            write_numeric_attribute(obj,dset_id,'size',sz);
            H5D.close(dset_id);
            H5S.close(space_id);
        end
        function write_numeric(obj,val,location)
            type = class(val); sz = size(val); count = prod(sz);
            if 0==count && verLessThan('matlab','8.5')
                space_id = H5S.create_simple(1,1,[]);
                dset_id = H5D.create(obj.fid,location,obj.(type),space_id,'H5P_DEFAULT');
                H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',0);
            else
                space_id = H5S.create_simple(ndims(val),fliplr(sz),[]);
                dset_id = H5D.create(obj.fid,location,obj.(type),space_id,'H5P_DEFAULT');
                H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',val);
            end
            write_char_attribute(obj,dset_id,'class',type);
            write_numeric_attribute(obj,dset_id,'size',sz);
            H5D.close(dset_id);
            H5S.close(space_id);
        end
        function val = read_variable(obj,location)
            id = H5O.open(obj.fid,location,'H5P_DEFAULT');
            info = H5O.get_info(id);
            type = read_char_attribute(obj,id,'class');
            sz = read_numeric_attribute(obj,id,'size');
            if strncmp(type,'ml',2), type = 'struct'; end
            switch info.type
                case 0
                    switch type
                        case 'struct'
                            nfield = read_numeric_attribute(obj,id,'nfield');
                            count = prod(sz) * (0~=nfield);
                            switch count
                                case 0, val = [];
                                case 1
                                    names = who(obj,location);
                                    for m=1:length(names), val.(names{m}) = read_variable(obj,[location '/' names{m}]); end
                                otherwise
                                    for m=1:count, val(m) = read_variable(obj,[location '/' num2str(m)]); end %#ok<AGROW>
                                    val = reshape(val,sz);
                            end
                        case 'cell'
                            count = prod(sz);
                            val = cell(sz);
                            if 0~=count
                                names = who(obj,location);
                                for m=1:length(names), val{m} = read_variable(obj,[location '/' names{m}]); end
                            end
                    end
                case 1
					if 0==prod(sz)  % If the size attribute is 0, do not read the dataset. See write_numeric().
						val = cast([],type);
					else
						val = H5D.read(id);
						switch type
							case 'char', val = val{1};
							case 'logical', val = logical(val);
						end
						val = reshape(val,sz);
					end
            end
            close(id);
        end
        function write_char_attribute(~,id,name,val)
            type_id = H5T.copy('H5T_C_S1');
            H5T.set_size(type_id,length(val));
            space_id = H5S.create('H5S_SCALAR');
            attr_id = H5A.create(id,name,type_id,space_id,'H5P_DEFAULT');
            H5A.write(attr_id,'H5ML_DEFAULT',val);
            H5A.close(attr_id);
            H5T.close(type_id);
        end
        function write_numeric_attribute(obj,id,name,val)
            type = class(val); sz = size(val);
            space_id = H5S.create_simple(ndims(val),fliplr(sz),fliplr(sz));
            attr_id = H5A.create(id,name,obj.(type),space_id,'H5P_DEFAULT');
            H5A.write(attr_id,'H5ML_DEFAULT',val);
            H5A.close(attr_id);
        end
        function val = read_numeric_attribute(~,id,name)
            attr_id = H5A.open(id,name);
            val = H5A.read(attr_id);
            H5A.close(attr_id);
        end
        function val = read_char_attribute(~,id,name)
            attr_id = H5A.open(id,name);
            val = H5A.read(attr_id)';
            H5A.close(attr_id);
        end
%         function write_cellstr(obj,val,location)
%             type_id = H5T.copy('H5T_C_S1');
%             H5T.set_size(type_id,'H5T_VARIABLE');
%             space_id = H5S.create_simple(1,length(val),[]);
%             dset_id = H5D.create(obj.fid,location,type_id,space_id,'H5P_DEFAULT');
%             H5D.write(dset_id,obj.string,'H5S_ALL','H5S_ALL','H5P_DEFAULT',val);
%             H5D.close(dset_id);
%             H5S.close(space_id);
%             H5T.close(type_id);
%         end
    end
end
