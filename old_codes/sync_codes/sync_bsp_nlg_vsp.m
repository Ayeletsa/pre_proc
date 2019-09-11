function p = sync_bsp_nlg_vsp(p)

bsp_dir =  fullfile(p.path_day_dir, 'bsp');
nlx_dir =  fullfile(p.path_day_dir, 'nlx');
sync_dir = fullfile(p.path_day_dir, 'sync');
out_dir = fullfile(sync_dir, 'sync_nlx_to_nlg');

% sync between NLG and BSP
time_conv_p_msec_nlx2bsp = sync_nlg_to_bsp(nlx_dir,bsp_dir,out_dir);
p.time_conv_p_msec_nlx2bsp = time_conv_p_msec_nlx2bsp;


%% if there are audio recordings, sync between VSP ts and BSP ts

if p.Audio == 1
    vsp_dir = fullfile(p.path_day_dir, 'vsp');
    out_dir = fullfile(sync_dir, 'sync_vsp_to_bsp');
    
    if ~exist(fullfile(out_dir,'vsp_ts_msec_fitted_to_bsp.mat'),'file') %if file does not exist
        % a. find the tone ts on the VSP recordings
        [vsp_ts_msecs, vsp_tone_ts_msecs] = Find_Tone_Times_in_Audio (vsp_dir,out_dir,date);
        % b. sync ts
        vsp_ts_msec_fitted_to_bsp = sync_vsp_to_bsp_through_nlg (nlx_dir,out_dir,vsp_tone_ts_msecs,vsp_ts_msecs,p.time_conv_p_msec_nlx2bsp);
        
    else %if file exist
        vsp_ts_msec_fitted_to_bsp = load(fullfile(out_dir,'vsp_ts_msec_fitted_to_bsp'));
        vsp_ts_msec_fitted_to_bsp = vsp_ts_msec_fitted_to_bsp.vsp_ts_msec_fitted_to_bsp;
    end
    
    p.vsp_ts_msec_fitted_to_bsp = vsp_ts_msec_fitted_to_bsp;
    
end

end