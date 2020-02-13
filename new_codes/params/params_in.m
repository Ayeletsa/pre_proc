%% Parameters to define for pre proc
bat='2336';
number_of_tt=16;
%% params file names:
p_in.spikes_params_file_name=['D:\Ayelet\Data\Data_Nlg_Proc\2bat_proj\spike_params_bat_',num2str(bat),'.mat'];
p_in.bsp_params_file_name=['D:\Ayelet\Data\Data_Nlg_Proc\2bat_proj\bsp_params_bat_',num2str(bat),'.mat'];
p_in.filter_params_file_name=['D:\Ayelet\Data\Data_Nlg_Proc\2bat_proj\filter_params_bat_',num2str(bat),'.mat'];

%%  1. excel file name:
%p_in.excel_sheet = 'D:\Matlab\pre_proc\test_pre_proc.xlsx';
p_in.excel_sheet = ['D:\Matlab\pre_proc\new_codes\inclusion_lists\BAT',num2str(bat),'_inclusion_list.xls'];

%% 2. excel rows for analysis
%p_in.day_rows = 1:9;
p_in.day_rows = 1:9;

p_in.cell_rows = 1:95;

%% 3. folder names:
p_in.path_datain = 'L:\Data\2batproj\yr_2019_bat_2336\'; % prefix for all data folders
p_in.path_dataout = 'L:\Data\2batproj\Data_Nlg_Proc\'; % prefix for output folders

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

p_in.spikes_params.thr =6; %check which threshold!
p_in.spikes_params.lib_corr_thr  = 0.8;
%p_in.spikes_params.coincidence_window=24;
p_in.spikes_params.CD_thr =6;
p_in.spikes_params.t_start_end = [];
p_in.spikes_params.thr_type = 'median'; % using median(abs(signal)) * factor
p_in.spikes_params.lib_spike_shapes = 'D:\Matlab\pre_proc\new_codes\neural_proc_codes\library_of_acceptable_spike_shapes_new.mat';
% params.CD_detect_win_len = 32;
p_in.spikes_params.CD_detect_win_len = 4;
p_in.spikes_params.CD_invalid_win_len = p_in.spikes_params.CD_detect_win_len*2;
% params.CD_n_ch_thr = 9;
% params.CD_n_TT_thr  = length(params.TT_to_use);
% params.CD_n_ch_thr = 0.5 * sum(params.active_TT_channels(:)); % at least on half of the channels
p_in.spikes_params.is_save_artifacts = 1;
p_in.spikes_params.number_of_wires_in_TT=4;
p_in.spikes_params.min_sep_events = 24;
p_in.spikes_params.nSamples = 32;
p_in.spikes_params.AlignSample = 8;
p_in.isolation_distance_in_min=10;

if number_of_tt==16;
    
    p_in.spikes_params.CD_n_TT_thr=6;
else
    p_in.spikes_params.CD_n_TT_thr=4;
end




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
p_in.filter_params.run_LFP=0;
p_in.filter_params.fwin = 2;         % we will run over the data in 2-min windows but save           % Same as we use for filtering ripples
p_in.filter_params.passband_spikes   = [600 6000];        % Filter for spikes
p_in.filter_params.passband_LFP      = [0.5 400];         % Filter for LFPs
p_in.filter_params.LFP_resamlpe_fs     = 2000;
p_in.filter_params.spikes_filter_type='highpassfir1';
% save  pramas
params_struct=p_in.filter_params;
params_file_name=p_in.filter_params_file_name;
save(num2str(params_file_name), '-struct', 'params_struct')

%% 7. bsp parameters

p_in.bsp_params.outliers_speed=20; %m/s
p_in.bsp_params.outliers_pos_dist_from_midline=2; %m
p_in.bsp_params.resample_fs = 100;
p_in.bsp_params.clib_file_name='D:\Matlab\pre_proc\new_codes\pos_codes\clib_files\20170809__tunnel_midline\calib_tunnel.mat'; %write the full path!\
p_in.bsp_params.landmark_file='D:\Matlab\pre_proc\new_codes\pos_codes\clib_files\Landmarks.xlsx';

% save  pramas
params_struct=p_in.bsp_params;
params_file_name=p_in.bsp_params_file_name;
save(num2str(params_file_name), '-struct', 'params_struct')

