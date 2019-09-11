function p=PRE_sync(p)
main_dir = p.path_day_dir;
sync_dir=fullfile(main_dir,'sync');

if ~exist(sync_dir)
    mkdir(sync_dir)
end
%% extract TTLs:
%a. bsp:
bsp_dir = fullfile(main_dir,'bsp');

bsp_TTL_ts_msec=initial_extract_bsp_TTL(bsp_dir);

%b. nlg
if p.nlg_self
    nlx_dir =  fullfile(p.path_day_dir, 'nlx');
    nlg_TTL_ts_msec=initial_extract_nlg_TTL(nlx_dir);
end

%b. aud
if p.Audio_self
    aud_self_dir =  fullfile(p.path_day_dir,p.Audio_dir_self);
    aud_self_TTL_ts_msec=initial_extract_nlg_TTL(aud_self_dir);
end
if p.Audio_other
    aud_other_dir =  fullfile(p.path_day_dir,p.Audio_dir_other);
    aud_other_TTL_ts_msec=initial_extract_nlg_TTL(aud_other_dir);
end

%% get sync ts:
%a. nlg
if p.nlg_self
    TTL_x=nlg_TTL_ts_msec;
    x_name='nlg';
    TTL_y=bsp_TTL_ts_msec;
    y_name='bsp';
    sync_name='sync_bsp_and_nlg';
    sync_ts=PRE_sync_x_2_y(TTL_x,TTL_y,x_name,y_name,sync_name,sync_dir);
    p.sync.bsp_ts_for_sync_with_nlg=sync_ts.Y;
    p.sync.nlg_ts_for_sync_with_bsp=sync_ts.X;
end

%b. aud
if p.Audio_self
    TTL_x=aud_self_TTL_ts_msec;
    x_name='aud_self';
    TTL_y=bsp_TTL_ts_msec;
    y_name='bsp';
    sync_name='sync_bsp_and_aud_self';
    sync_ts=PRE_sync_x_2_y(TTL_x,TTL_y,x_name,y_name,sync_name,sync_dir);
    p.sync.bsp_ts_for_sync_with_aud_self=sync_ts.Y;
    p.sync.aud_self_ts_for_sync_with_bsp=sync_ts.X;
end


if p.Audio_other
    TTL_x=aud_other_TTL_ts_msec;
    x_name='aud_other';
    TTL_y=bsp_TTL_ts_msec;
    y_name='bsp';
    sync_name='sync_bsp_and_aud_other';
    sync_ts=PRE_sync_x_2_y(TTL_x,TTL_y,x_name,y_name,sync_name,sync_dir);
    p.sync.bsp_ts_for_sync_with_aud_other=sync_ts.Y;
    p.sync.aud_other_ts_for_sync_with_bsp=sync_ts.X;
end

%% sync nlg events if needed:
if strcmp(p.sync_to,'bsp') && ~isempty(p.S)
    p.S.start_time=interp1(p.sync.nlg_ts_for_sync_with_bsp,p.sync.bsp_ts_for_sync_with_nlg, p.S.start_time*1e3, 'linear','extrap');
    p.S.end_time=interp1(p.sync.nlg_ts_for_sync_with_bsp,p.sync.bsp_ts_for_sync_with_nlg, p.S.start_end*1e3, 'linear','extrap');

end
