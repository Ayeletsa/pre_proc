function PRE_create_day_structs (p_in)
% Function description:
% Creates day struct with all the neural data and bsp data.
% If there is already a day struct for this day the function will continue to the next day
% For days with audio recordings it can create also audio struct data


excel_sheet=p_in.excel_sheet;
day_rows=p_in.day_rows;
%% 1. initial extract:
% a. Transform from nlg files to nlx files (neuro and audio)
PRE_loop_nlg2nlx (excel_sheet,'Experiments',day_rows)

%% 2. Read excel table for desired recording days

P = NLG_PRE_read_and_fill_excel_sheet (p_in,excel_sheet,'Experiments',day_rows);
day_struct_dir=fullfile(P(1).path_dataout,P(1).year_bat_path,'day_structs');

if ~exist(day_struct_dir)
    mkdir(day_struct_dir)
end

for ii_rec = 1 : length(day_rows)
    close all
    p = P(ii_rec);
    struct_name=fullfile(p.path_dataout,p.year_bat_path,'day_structs',['bat_',num2str(p.bat),'_day_',num2str(p.day),'.mat']);
    
    if ~exist(struct_name)
        %3. sync:
        p=PRE_sync(p);
        
        if p.nlg_self
        % 3. filter for spikes and LFP:
        % ------------------------------
        PRE_filter_CSCs(p) %note this is not synced!
       
        % 4. detect spikes:
        % ------------------------------
        p=Nlx_detect_spikes_CSC3(p);
        %p = Nlg_detect_spikes_library_first_findpeaks_1(p);
        end
        % 5. extract BSP data:
        % ------------------------------

        p = BSP_PRE_PROC_data (p); %TO DO -  take parameters to the param_in script
         
      
        % 6. save day struct:
        % ------------------------------
        save(struct_name,'p')
   end 
        if p.Audio_self == 1 || p.Audio_other == 1%only if there are audio recordings:
            load(struct_name)
        % 7. create audio strcut:
        % ------------------------------        
            PRE_create_audio_struct(p)
        end
    
end


end

