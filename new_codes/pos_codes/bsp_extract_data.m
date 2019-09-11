function [bsp_data] = bsp_extract_data(bsp_dir)
%%% written by Tamir, changed by Shir 04jan2017

%% read files
dist_csv_file_names = dir(fullfile(bsp_dir,['*_distances.csv']));
TTL_ts_ns_all = [];
bsp_data = struct('tag_ID',[], 'pos', [], 'ts_ns', []);
for ii_file=1:length(dist_csv_file_names)
    file_name = fullfile(bsp_dir,dist_csv_file_names(ii_file).name);
    ii_file
    % replace ';' with ','
    fid  = fopen(file_name,'r');
    f=fread(fid,'*char')';
    fclose(fid);
    f = strrep(f,';',',');
    file_name = strrep(file_name, 'distances.csv', 'distances_COMMA.csv');
    fid  = fopen(file_name,'w');
    fprintf(fid,'%s',f);
    fclose(fid);
    
    [NUM,TXT,RAW]=xlsread(file_name);
    
    %% read header
    num_tags = RAW{2,2};
    for ii_tag=1:num_tags
        clear tag_id ts_ns X Y Z pos
        tag_id = RAW{3,(ii_tag-1)*3+2};
        tag_id = tag_id(6:end);
        tag_id = str2num(tag_id);
        ts_ns = NUM(4:end,1);
        X = NUM(4:end,(ii_tag-1)*3+2);
        Y = NUM(4:end,(ii_tag-1)*3+3);
        Z = NUM(4:end,(ii_tag-1)*3+4);
        pos = [X Y Z];
        nan_IX = find(isnan(X));
        pos(nan_IX,:) = [];
        ts_ns(nan_IX) = [];
        
        % insert data to relevant tag (if tag entry does not exist, create new one)
        if ismember(tag_id, [bsp_data.tag_ID])
            tag_IX = find([bsp_data.tag_ID] == tag_id);
        else
            tag_IX = length([bsp_data.tag_ID])+1;
            bsp_data(tag_IX).tag_ID = tag_id;
        end
        bsp_data(tag_IX).pos = [bsp_data(tag_IX).pos; pos];
        bsp_data(tag_IX).ts_ns = [bsp_data(tag_IX).ts_ns; ts_ns];
    end
    
    % add TTL (if exist)
    TTL_col = 1+num_tags*3+1;
    if size(NUM,2)>= TTL_col
        TTL_IX = find(NUM(4:end, TTL_col) == 1);
        ts_ns = NUM(4:end,1);
        TTL_ts_ns = ts_ns(TTL_IX);
        TTL_ts_ns_all = [TTL_ts_ns_all; TTL_ts_ns];
    else
        warning('No TTL data in current file, CHECK WHY!!!!')
    end
    
end

%% save as mat file
bsp_TTL_ts_ns = TTL_ts_ns_all;
save(fullfile(bsp_dir,'bsp_data'), 'bsp_data');
save(fullfile(bsp_dir,'bsp_TTL'), 'bsp_TTL_ts_ns');
for ii_tag = 1:length([bsp_data.tag_ID])
    tag_id = bsp_data(ii_tag).tag_ID;
    bsp_pos = bsp_data(ii_tag);
    save(fullfile(bsp_dir,['bsp_pos_tag_' num2str(tag_id)]), 'bsp_pos');
end

end