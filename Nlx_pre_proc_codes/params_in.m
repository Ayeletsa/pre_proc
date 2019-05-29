%% Parameters to define for pre proc
bat='2389';
%% params file names:
p_in.spikes_params_file_name=['D:\Ayelet\Data\Data_Nlg_Proc\2bat_proj\spike_params_bat_',num2str(bat),'.mat'];
p_in.bsp_params_file_name=['D:\Ayelet\Data\Data_Nlg_Proc\2bat_proj\bsp_params_bat_',num2str(bat),'.mat'];
p_in.filter_params_file_name=['D:\Ayelet\Data\Data_Nlg_Proc\2bat_proj\filter_params_bat_',num2str(bat),'.mat'];

%%  1. excel file name:
p_in.excel_sheet = 'D:\Matlab\pre_proc\test_pre_proc.xlsx';
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

%tamir's params:
p_in.spikes_params.use_neg_thr =0; %find also negative events
p_in.spikes_params.thr =repmat(7,4,4); %check which threshold!
p_in.spikes_params.lib_corr_thr  = 0.8;
p_in.spikes_params.coincidence_window=24;
p_in.spikes_params.CD_thr =repmat(6,4,4);
p_in.spikes_params.t_start_end = [];
p_in.spikes_params.thr_type = 'median'; % using median(abs(signal)) * factor
p_in.spikes_params.lib_spike_shapes = 'D:\Matlab\pre_proc\library_of_acceptable_spike_shapes_new.mat';
% params.CD_detect_win_len = 32;
p_in.spikes_params.CD_detect_win_len = 4;
p_in.spikes_params.CD_invalid_win_len = 32*2;
% params.CD_n_ch_thr = 9;
% params.CD_n_TT_thr  = length(params.TT_to_use);
% params.CD_n_ch_thr = 0.5 * sum(params.active_TT_channels(:)); % at least on half of the channels
p_in.spikes_params.is_save_artifacts = 1;
p_in.spikes_params.number_of_wires_in_TT=4;
p_in.spikes_params.min_sep_events = 24;
p_in.spikes_params.nSamples = 32;
p_in.spikes_params.AlignSample = 8;

% save spikes pramas
params_struct=p_in.spikes_params;
params_file_name=p_in.spikes_params_file_name;
save(num2str(params_file_name), '-struct', 'params_struct')



% %old params:
% p_in.library_file_name='D:\Matlab\Spike_detection_proj\Data_for_library\new_lib.mat';
% p_in.include_negative_threshold=1; %find also negative events
% p_in.do_coincidence_detection=1;
% p_in.compare_to_library=1;
% 
% p_in.threhold_factor=5; %for voltage threshold (threhold=median*factor)
% p_in.min_sep_events=15; %min separation between events - ~=500us
% p_in.r_threshold = 0.9;
% p_in.coincidence_window=500; % Duration of coincidence-detection window (in microsec) ACROSS tetrodes
% p_in.thresh_for_coinc=4; %remove if X tetrodes has event in the window
% 
% p_in.create_movie=0;

%% 6. filtering parameters

%% 7. bsp parameters

p_in.bsp_params.outliers_speed=20; %m/s
p_in.bsp_params.outliers_pos_dist_from_midline=2; %m
p_in.bsp_params.resample_fs = 100;
p_in.bsp_params.clib_file_name='D:\Ayelet\Data\2batproj\20170809__tunnel_midline\calib_tunnel.mat'; %write the full path!\
p_in.bsp_params.landmark_file='D:\Ayelet\Data\2batproj\Landmarks.xlsx';

% save spikes pramas
params_struct=p_in.bsp_params;
params_file_name=p_in.bsp_params_file_name;
save(num2str(params_file_name), '-struct', 'params_struct')

