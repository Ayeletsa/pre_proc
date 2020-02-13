function PP = NLG_PRE_read_and_fill_excel_sheet (p_in,excel_sheet,sheet,rows)

%% Read excel sheet and create parameter structure from it

P = PRE_read_excel_sheet(excel_sheet,sheet,rows,p_in.numeric_fields,p_in);


if length(P) == 0 %if no records are entered
    return
end


%% Extract events for each day
records = 1:length(P);
% records(ismember(records,[1,12,15,16,19,21,24]))=[];
for ii_rec = records % number of records (days) in excel
    
    p = P(ii_rec);
    
    %Generating path for the data folder
    day_string= num2str(p.day);
    p.path_bat_date=['BAT' num2str(p.bat) '\' day_string '\nlg\'];
    
    disp('========================================================');
    disp(' ');
    disp(['Extracting day' num2str(ii_rec) ' ' day_string]);
    disp(' ');
    
    
    %----------------------------------------------------------------------
    % 1) Extracting ans Show Events and Timestamps for both Nlg and Nlx data
    %----------------------------------------------------------------------
    
    
    %Neurologger: Extracting Events, Event Timestamps, and TTL times from the event log of the Nlg
   if p.nlg_self
       if p.take_event_from_audio==1
           file_name_Event_Nlg=  fullfile(p.path_day_dir,p.Audio_dir_other,'\');
           
       else
           
           file_name_Event_Nlg=  [ p.path_day_dir '\nlx\'];
       end
    [p.Nlg_EventStrings, p.Nlg_EventTimestamps] = NLG_PRE_read_excel_events_and_TTL(file_name_Event_Nlg);
 
   
    disp('Neurologger EVENT LIST:');
    disp('===========');
    p.Nlg_EventStrings % displaying event list
    
    
    %----------------------------------------------------------------------
    % 2) Check that all info is filled in in excel
    %----------------------------------------------------------------------
    
    %1) check that all Nlg sessions are defined in excel
    if ~isfield(p,'S') || isempty(p.S)
        disp('please complete Nlg session names info in excel and run again');
        return;
    end
    %2) check that all Nlg sessions events are defined in excel
    nsessions = length(p.S);
    for nses = 1:nsessions
        s = p.S(nses);
        if isempty(s.events) || all(isnan(s.events))
            disp('please complete Nlg session events info in excel and run again');
            return;
        end
    end
    
    
    %----------------------------------------------------------------------
    % 3) For each session defined in excel get the following parameters:
    %----------------------------------------------------------------------
    
    % start_behavior  - pointer to start of behavior in event list
    % end_behavior    - pointer to end of behavior in event list
    % event_list      - event list itself
    % time_stamps     - time stamps of each event in event list
    % time_stamps offsets - if entered manually by the user
    
    % Important:
    %Nlx - timestamps in microsec
    %Nlg - timestamps also in microseconds (!)
    
    p = NLG_PRE_get_session_times_Nlg(p);
    
   end 
    %----------------------------------------------------------------------
    % 9) Update the structure (inclusion list)
    %----------------------------------------------------------------------
    PP(ii_rec) = p;
    
end % loops on excel records for display of events  only



end