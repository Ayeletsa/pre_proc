function p = PRE_sync_aud2nlx(p,aud_dir,self)

aud_dir =  fullfile(p.path_day_dir, aud_dir);
nlx_dir =  fullfile(p.path_day_dir, 'nlx');
sync_dir = fullfile(p.path_day_dir, 'sync');
out_dir = fullfile(sync_dir, 'sync_aud2nlx');

if ~exist(out_dir,'dir')
    mkdir(out_dir)
end

%% aud TTL
aud_TTL_file_name = fullfile(aud_dir, 'EVENTS__Digital in.nev');
FieldSelection = [1 0 0 0 0];
ExtractHeader = 0;
ExtractMode = 1;
ModeArray = [];
aud_TTL_ts_usec = Nlx2MatEV( aud_TTL_file_name ,FieldSelection,ExtractHeader,ExtractMode,ModeArray);
aud_TTL_ts_msec = aud_TTL_ts_usec*1e-3;
aud_TTL_intervals = diff(aud_TTL_ts_msec);
aud_TTL_intervals_inc = diff(aud_TTL_intervals);

%% NLG TTL
nlg_TTL_file_name = fullfile(nlx_dir, 'EVENTS__Digital in.nev');
FieldSelection = [1 0 0 0 0];
ExtractHeader = 0;
ExtractMode = 1;
ModeArray = [];
nlg_TTL_ts_usec = Nlx2MatEV( nlg_TTL_file_name ,FieldSelection,ExtractHeader,ExtractMode,ModeArray);
nlg_TTL_ts_msec = nlg_TTL_ts_usec*1e-3;
nlg_TTL_intervals = diff(nlg_TTL_ts_msec);
nlg_TTL_intervals_inc = diff(nlg_TTL_intervals);

%% plot TTL timings
ax_h = [];
figure('units','normalized','outerposition',[0 0 1 1])
ax_h(1) = subplot(2,2,1);
plot(aud_TTL_ts_msec(2:end)*1e-3/60,aud_TTL_intervals*1e-3,'.')
% xlabel('Time (minutes)')
ylabel('inter-TTL-interval (sec)')
title('aud TTL timings')
ax_h(3) = subplot(2,2,3);
plot(aud_TTL_ts_msec(3:end)*1e-3/60,aud_TTL_intervals_inc,'.')
xlabel('Time (minutes)')
ylabel('inter-TTL-interval increament (msec)')
ax_h(2) = subplot(2,2,2);
plot(nlg_TTL_ts_msec(2:end)*1e-3/60,nlg_TTL_intervals*1e-3,'.')
title('nlg TTL timings')
% xlabel('Time (minutes)')
ylabel('inter-TTL-interval (sec)')
ax_h(4) = subplot(2,2,4);
plot(nlg_TTL_ts_msec(3:end)*1e-3/60,nlg_TTL_intervals_inc,'.')
xlabel('Time (minutes)')
ylabel('inter-TTL-interval increament (msec)')
ylimits = get(ax_h([1 2]), 'ylim');
ylimits = [  min([ylimits{:}]) max([ylimits{:}])  ];
linkaxes(ax_h([1 2]),'y')
set(ax_h(1), 'ylim' , ylimits);
ylimits = get(ax_h([3 4]), 'ylim');
ylimits = [  min([ylimits{:}]) max([ylimits{:}])  ];
linkaxes(ax_h([3 4]),'y')
set(ax_h(3), 'ylim' , ylimits);
saveas(gcf, fullfile(out_dir, 'sync_aud_nlg__TTL_timinigs'), 'jpeg')
saveas(gcf, fullfile(out_dir, 'sync_aud_nlg__TTL_timinigs'), 'fig')
close(gcf)

%% sync
N_smoothing=2; %check these parameters!
max_diff=100; %check these parameters!
time_conv_p_msec_aud2nlx = sync_TTL_polyfit(round(aud_TTL_ts_msec),round(nlg_TTL_ts_msec), N_smoothing, max_diff); % p values to convert nlg ts to bsp tssaveas(gcf, fullfile(out_dir, 'sync_bsp_nlg'), 'jpeg')
close(gcf)
save( fullfile(out_dir, 'time_conv_p_msec_bsp2nlx') , 'time_conv_p_msec_bsp2nlx')
if self==1
    p.time_conv_p_msec_aud2nlx_self = time_conv_p_msec_aud2nlx;
else
    p.time_conv_p_msec_aud2nlx_other = time_conv_p_msec_aud2nlx;
end
end
