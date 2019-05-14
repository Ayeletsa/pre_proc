function [vsp_ts_msec_fitted_to_bsp,pp,chosen_slope] = sync_vsp_to_bsp_through_nlg (main_dir,vsp_tone_ts_msecs,vsp_time_vec_msecs,time_conv_p_msec_nlx2bsp)

nlx_dir =  fullfile(main_dir, 'nlx');
sync_dir = fullfile(main_dir, 'sync');
out_dir = fullfile(sync_dir, 'sync_vsp2bsp');

%% read tone timestamps 

nlx_tones_file_name = fullfile(nlx_dir, 'EVENTS__Tone generated.nev');
FieldSelection = [1 0 0 0 0];
ExtractHeader = 0;
ExtractMode = 1;
ModeArray = [];
nlx_tone_ts_usec = Nlx2MatEV( nlx_tones_file_name ,FieldSelection,ExtractHeader,ExtractMode,ModeArray);
nlx_tone_ts_msecs = nlx_tone_ts_usec* 1e-3;

%% sync VSP to NLG

N_smoothing = 15;
max_diff = 15000;

[vsp_ts_msec_fitted_to_nlg,pp,chosen_slope] = sync_nlg_to_vsp_with_inital_fit(vsp_tone_ts_msecs, nlx_tone_ts_msecs, N_smoothing, max_diff, vsp_time_vec_msecs, out_dir);

%% sync VSP to BSP

vsp_ts_msec_fitted_to_bsp = polyval(time_conv_p_msec_nlx2bsp,vsp_ts_msec_fitted_to_nlg);
save(fullfile(out_dir,'vsp_ts_msec_fitted_to_bsp'),'vsp_ts_msec_fitted_to_bsp','-v7.3')

end
