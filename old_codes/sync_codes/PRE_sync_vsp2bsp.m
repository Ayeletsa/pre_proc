function p = PRE_sync_vsp2bsp (p,sigAudioUnfilt)

vsp_dir = fullfile(p.path_day_dir, 'vsp');
nlx_dir =  fullfile(p.path_day_dir, 'nlx');
sync_dir = fullfile(p.path_day_dir, 'sync');
out_dir = fullfile(sync_dir, 'sync_vsp2bsp');


%% if there are audio recordings, sync between VSP ts and BSP ts

    % a. find the tone ts on the VSP recordings
    [vsp_ts_msecs, vsp_tone_ts_msecs] = Find_Tone_Times_in_Audio (sigAudioUnfilt,out_dir);

    % b. sync ts
    [vsp_ts_msec_fitted_to_bsp,time_conv_pp_msec_vsp2nlg,slope_for_vsp_ts] = sync_vsp_to_bsp_through_nlg (nlx_dir,out_dir,vsp_tone_ts_msecs,vsp_ts_msecs,p.time_conv_p_msec_nlx2bsp);
    
end
