%% Parameters to define for pre proc
%%  1. excel file name:
p_in.excel_sheet = 'D:\Matlab\pre_proc\BAT2389_inclusion_list.xls';
%p_in.excel_sheet = 'D:\Matlab\pre_proc\BAT2287_inclusion_list.xls';

%% 2. excel rows for analysis
p_in.day_rows = 1:18;
p_in.cell_rows = 1:119;

%% 3. folder names:
p_in.path_datain = 'D:\Ayelet\Data\2batproj\'; % prefix for all data folders
p_in.path_dataout = 'D:\Ayelet\Data\Data_Nlg_Proc\'; % prefix for output folders

%% 4. excel file fields:
p_in.numeric_fields = {'bat','day','reference_channel','TT','depth', ...
    'active_channels','nsessions','use_for_sorting', 'use_for_analysis', ...
    'use_tetrodes','throw_away_times', ...
    'events_#','time_offsets_#_in_seconds','exact_time_in_microseconds_#' 'Nlx_events', 'Nlx_time_offsets_in_seconds','Nlx_exact_time_in_microseconds','Sync_by_LED',...
    'cell_num', 'cell_id','b1_similar_cell','b3_similar_cell', 'single_unit','pyramidal', 'behaviorally_active'}; % numeric fields in excel
p_in.Nlg_EventStrings=[]; % will be filled automatically later on
p_in.Nlg_EventTimestamps=[]; % will be filled automatically later on

%% 5. spike detection parametes:
p_in.library_file_name='D:\Matlab\Spike_detection_proj\Data_for_library\new_lib.mat';

p_in.include_negative_threshold=1; %find also negative events
p_in.do_coincidence_detection=1;
p_in.compare_to_library=1;

p_in.threhold_factor=5; %for voltage threshold (threhold=median*factor)
p_in.min_sep_events=15; %min separation between events - ~=500us
p_in.r_threshold = 0.9;
p_in.coincidence_window=500; % Duration of coincidence-detection window (in microsec) ACROSS tetrodes
p_in.thresh_for_coinc=4; %remove if X tetrodes has event in the window

p_in.create_movie=0;
