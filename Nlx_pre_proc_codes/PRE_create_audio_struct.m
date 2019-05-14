function PRE_create_audio_struct(p)

struct_folder = fullfile(p.path_dataout,p.year_bat_path,'day_structs');
struct_name = sprintf('bat_%d_day_%d_audio_struct',p.bat,p.day);
file_name = fullfile(struct_folder,struct_name);

if ~exist(file_name,'file')
    %1. convert to nlx
    self_dir=fullfile(p.path_day_dir,p.Audio_dir_self) ;
    nlg2nlx_audio(self_dir)
    if p.Audio_dir_other~=0
        other_dir=fullfile(p.path_day_dir,p.Audio_dir_other) ;
        nlg2nlx_audio(other_dir)
    end
    %2. read the data
    Filename=fullfile(self_dir,'audio.ncs');
    [signal, ts, fs] = Nlx_csc_read(Filename,[]);
    p.Aud.self_signal=signal;
    p.Aud.self_ts=ts;
    p.Aud.fs=fs;
    if p.Audio_dir_other~=0
        
        Filename=fullfile(other_dir,'audio.ncs');
        [signal, ts, fs] = Nlx_csc_read(Filename,[]);
        p.Aud.other_signal=signal;
        p.Aud.other_ts=ts;
    end
    %3. sync to nlg
    self=1;
    p=PRE_sync_aud2nlx(p,self_dir,self);
    
    if p.Audio_dir_other~=0
        self=0;
        p=PRE_sync_aud2nlx(p,other_dir,self);
    end
    A = p;
    A.sigAudioU = signal;
    A.vsp_ts_msec_fitted_to_bsp = vsp_ts_msec_fitted_to_bsp;
    A.time_conv_pp_msec_vsp2nlg = time_conv_pp_msec_vsp2nlg;
    A.slope_for_vsp_ts = slope_for_vsp_ts;
    
    save(file_name,'A','-v7.3')
    
end

end