function p = NLG_PRE_loop_Nlg2spikes (p)

main_dir = p.path_day_dir;
active_channels = p.active_channels;
PRE_filter_CSCs(main_dir,active_channels)
active_TTs = p.use_tetrodes;
p = Nlg_detect_spikes_new_code (p,main_dir,active_TTs,active_channels);

end