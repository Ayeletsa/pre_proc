function [event_time_us,TTL_string] = get_events_from_local_log(main_dir,text)
txt_files = dir([main_dir '\*.txt']);

TTL_string = {};
event_time_us = [];
for i_file = 1:numel(txt_files)
    local_log = txt_files(i_file).name;
    fid = fopen([main_dir filesep local_log]);
    data=textscan(fid,'%s %f %s','Delimiter',',');
    fclose(fid);
    prefix_TTL = 't=';

     
    ttl_ind = cellfun(@(x) contains(x,text),data{3});
    TTL_string = [TTL_string;data{3}(ttl_ind)];
       TTL_pos = cellfun(@(x) strfind(x,prefix_TTL),TTL_string);

    TTL_time_us = cellfun(@(x,y) milliseconds(duration(x(y + length(prefix_TTL) : y + length(prefix_TTL) + 15))),TTL_string,num2cell(TTL_pos))*1e3;
    

   
    event_time_us = [event_time_us;data{2}(ttl_ind) * 1e3];
end

[event_time_us,ia,~] = unique(event_time_us);
event_time_us = event_time_us';
event_string = event_string(ia)';
end