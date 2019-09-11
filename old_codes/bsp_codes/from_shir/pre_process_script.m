% preprocessing script:
clear all
bat = '9845';
dataPath = ['D:\Project\DATA\' bat '\'];
files = [20170606];%[20170920 20170924 20170925 20170926 20170927 20170928 20171001 ...
%20171002 20171003 20171005 20171008 20171009 20171012 20171016];
tunnel_calib_path = 'D:\Project\Tamir\yr2017_bat0148_behav_neural\calibrations\20170809__tunnel_midline\calib_tunnel.mat';

%% Nlg2Nlx:
% for i=1:length(files)
%     files(i)
%     path = [dataPath num2str(files(i)) '\'];
%     Nlg2Nlx(path);
%     eval(['e' num2str(files(i)) '_bat' bat]);
%     Nlg2Nlx(path,param);
%     
%     clear param
%     close all
% end
%% extract BSP data:
% for i=1:length(files)
%     main_dir = [dataPath num2str(files(i)) '\bsp\'];
%     bsp_extract_data(main_dir);
% end
% %% sync BSP DATA to neurologger: 
% for i=1:length(files)
%     bsp_dir = [dataPath num2str(files(i)) '\bsp\'];
%     nlg_dir = [dataPath num2str(files(i)) '\nlx\'];
%     PRE_sync_bsp_to_nlg(bsp_dir,nlg_dir,bsp_dir);
% end
%% Linearize BSP data (including the removal of outlier position points):
for i = 1:length(files)
    expname = ['e' num2str(files(i)) '_bat' bat];
    eval(expname);
    load([param.path.bsp '\bsp_pos_tag_' num2str(param.general.BspTag) '.mat']);
    %%% for debugging:
    bsp_pos.ts_nlg_usec = bsp_pos.ts_ns*10^-3;
    %%%
    bsp_pos = linearize_bsp(bsp_pos,tunnel_calib_path);
    
    %save([param.path.bsp '\bsp_pos_tag_' num2str(param.general.BspTag) '_linearized.mat'],bsp_pos);
end
%% fill BSP data holes:
for i = 1:length(files)
    expname = ['e' num2str(files(i)) '_bat' bat];
    eval(expname);
    %load([param.path.bsp '\bsp_pos_tag_' num2str(param.general.BspTag) '_linearized.mat']);
    raw_data_pos = bsp_pos.pos_linear;
    raw_data_ts_usec = bsp_pos.ts_nlg_usec;
    [bsp_full_data_pos,bsp_full_data_ts] = fill_bsp_holes_final(raw_data_pos',raw_data_ts_usec');
    bsp_pos.pos_linear = bsp_full_data_pos;
    bsp_pos.ts_nlg_usec = bsp_full_data_ts;
    %save([param.path.bsp '\bsp_pos_tag_' num2str(param.general.BspTag) '_linearized.mat'],bsp_pos);
end
%% upsample BSP data:
for i = 1:length(files)
    expname = ['e' num2str(files(i)) '_bat' bat];
    eval(expname);
    %load([param.path.bsp '\bsp_pos_tag_' num2str(param.general.BspTag) '_linearized.mat']);
    bsp_pos_new = upsample_100hz_final(bsp_pos);
    %save([param.path.bsp '\bsp_pos_tag_' num2str(param.general.BspTag) '_linearized.mat'],bsp_pos);
end

%% filtering and spikes detection:
% for i=1:length(files)
%     expname = ['e' num2str(files(i)) '_bat' bat];
%     eval(expname);
%     disp(['expname: ' num2str(files(i))]);
%     Nlx_filter_CSCs(expname,1);
%     Nlg_detect_spikes_library_first_findpeaks(param);
% end
