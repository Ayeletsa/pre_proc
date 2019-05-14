function PRE_run_exp(exp_ID)

%%
exp_data_file = 'L:\Analysis\Code\inclusion_lists\recordings_summary.xlsx';
exp_t = readtable(exp_data_file, 'Sheet', 'Experiments', 'ReadRowNames',1);
exp = exp_t{exp_ID,:};

%% generaet data path strings
date_str = datestr(datevec(exp.date{1},'dd/mm/yyyy'), 'yyyymmdd' );
batNum_str = sprintf('%04d',exp.batNum);
clear exp_path
exp_path.bsp = fullfile('L:\DATA', [batNum_str '_' exp.batName{1}], date_str, 'bsp', 'client');
exp_path.nlg = fullfile('L:\DATA', [batNum_str '_' exp.batName{1}], date_str, 'nlg');
exp_path.nlx = fullfile('L:\DATA', [batNum_str '_' exp.batName{1}], date_str, 'nlx');
exp_path.nlx = fullfile('L:\DATA', [batNum_str '_' exp.batName{1}], date_str, 'nlx');
exp_path

%%
addpath(genpath('L:\Analysis\Code'))
% expname = 'exp_b0079_d20160919';
% addpath(genpath('D:\Tamir\PROJECTS\Neurologger\testing\Analysis\'))
% expname = 'mouse_logger_20170102';
% eval(expname);

%% Position (bsp) related
% bsp_extract_data(param.path.bsp)
% bsp_pre_process(param.path.bsp) % TODO: create
% bsp_calib_

%% Neural (Spikes+LFP)
% Nlg2Nlx(param.path.raw) % TODO: insert params from here, rather than change them inside Nlg2Nlx function
% Nlx_filter_CSCs(expname,0) % TODO: create a separate function PRE_filter_CSCs that will call Nlx_filter_CSCs (Nlx_filter_CSCs should be independant)
% Nlx_detect_spikes_CSC(expname,0)
% TODO: add some analysis AFTER spike sorting (separation matrices/FR/AC/xcorr/ISI/...)

%% sync position/neural
% PRE_sync_bsp_to_nlg(expname)
% TODO: make sure to plot figures showing good sync!
