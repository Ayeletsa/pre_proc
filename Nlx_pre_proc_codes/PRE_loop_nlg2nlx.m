function PRE_loop_nlg2nlx (excel_sheet,sheet,rows)
%Reads the excel and calls the function bellow for all relevant days in a loop. 
%It is important to check the ref channel that is updated in the excel


[num,txt,raw] = xlsread(excel_sheet,sheet);

field_names = txt(1,:);
nfields = length(field_names);
ref_channel_ind = find(strcmp(field_names,'reference_channel'));
path_day_dir_ind = find(strcmp(field_names,'path_day_dir'));
active_channels_ind = find(strcmp(field_names,'active_channels'));
active_TTs_ind = find(strcmp(field_names,'use_tetrodes'));

sheet = raw(2:end,:);
nrows = size(sheet,1);
rows(rows > nrows) = [];
sheet = sheet(rows,1:nfields);

for ii_rec = 1:size(sheet,1)

    ref_channel = cell2mat(sheet(ii_rec,ref_channel_ind));

    if ~isnumeric(ref_channel)
        disp ('Fill refernce channel correctly')
        return
    end

    main_dir = char(sheet(ii_rec,path_day_dir_ind));
    if ~ischar(main_dir)
        disp ('Fill path correctly')
        return
    end

Nlg2Nlx_tamir(main_dir)
ii_rec
end

end