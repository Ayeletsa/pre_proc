function PRE_create_day_structs (p_in)

excel_sheet=p_in.excel_sheet;
day_rows=p_in.day_rows;
%% transform raw recordings into filterd LFP and detected spikes
PRE_loop_nlg2nlx (excel_sheet,'Experiments',day_rows)

%% Read excel table for desired recording days

P = NLG_PRE_read_and_fill_excel_sheet (p_in,excel_sheet,'Experiments',day_rows);
day_struct_dir=fullfile(P(1).path_dataout,P(1).year_bat_path,'day_structs');

if ~exist(day_struct_dir)
    mkdir(day_struct_dir)
end

for ii_rec = 1 : length(day_rows)
    
    p = P(ii_rec);
    struct_name=fullfile(p.path_dataout,p.year_bat_path,'day_structs',['bat_',num2str(p.bat),'_day_',num2str(p.day),'.mat']);
    
    if ~exist(struct_name)
        PRE_filter_CSCs(p)
        p = Nlg_detect_spikes_library_first_findpeaks_1(p);
        p = BSP_PRE_extract_data (p); %TO DO - add landmarks  + take parameters to the param_in script
        p = PRE_sync_nlg2bsp(p);
        
        save(struct_name,'p')
    
        if p.Audio == 1
           
            PRE_create_audio_struct(p)
        end
    end
end


end

