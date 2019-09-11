function PRE_create_audio_struct(p)

struct_folder = fullfile(p.path_dataout,p.year_bat_path,'day_structs');
struct_name = sprintf('bat_%d_day_%d_audio_struct',p.bat,p.day);
file_name = fullfile(struct_folder,[struct_name,'.mat']);
self_dir=fullfile(p.path_day_dir,p.Audio_dir_self) ;

if ~exist(file_name,'file')
    %1. convert to nlx the other bat
    if p.Audio_dir_other~=0
        other_dir=fullfile(p.path_day_dir,p.Audio_dir_other) ;
        nlg2nlx_audio(other_dir)
    end
    %2. read the data
    Filename=fullfile(self_dir,'audio.ncs');
    [signal, ts, fs] = Nlx_csc_read(Filename,[]);
    p.Aud.self_signal=signal;
    % sync:
    if strcmp(p.sync_to,'bsp')
        ts=interp1(p.sync.aud_self_ts_for_sync_with_bsp,p.sync.bsp_ts_for_sync_with_aud_self, ts*1e3, 'linear','extrap');
    end
    p.Aud.self_ts=ts;
    p.Aud.fs=fs;
    if p.Audio_dir_other~=0
        Filename=fullfile(other_dir,'audio.ncs');
        [signal, ts, fs] = Nlx_csc_read(Filename,[]);
        if strcmp(p.sync_to,'bsp')
            ts=interp1(p.sync.aud_other_ts_for_sync_with_bsp,p.sync.bsp_ts_for_sync_with_aud_other, ts*1e3, 'linear','extrap');
        end
        p.Aud.other_signal=signal;
        p.Aud.other_ts=ts;
    end
    
%     % 3. sync to bsp
%     p.Aud.self_ts=PRE_sync_aud_2_bsp(p,self_dir,p.Aud.self_ts)/1e3;
%     p.Aud.other_ts=PRE_sync_aud_2_bsp(p,other_dir,p.Aud.other_ts)/1e3;
%     
    
    
    %     %3. sync to nlg
    %     if p.nlg
    %     self=1;
    %     p=PRE_sync_aud2nlx(p,self_dir,self);
    %
    %     if p.Audio_dir_other~=0
    %         self=0;
    %         p=PRE_sync_aud2nlx(p,other_dir,self);
    %     end
    %     end
    %
    
    
    save(file_name,'p','-v7.3')
    
end

end