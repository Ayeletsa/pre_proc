function   ts_synced=PRE_sync_aud_2_bsp(p,aud_dir,ts_to_convert)

bsp_dir =  fullfile(p.path_day_dir, 'bsp');
sync_dir = fullfile(p.path_day_dir, 'sync');
out_dir = fullfile(sync_dir, 'sync_aud2bsp');

if ~exist(out_dir,'dir')
    mkdir(out_dir)
end

% TODO add this option if there will be bugs!
sync_jump_ts=[];
 


%% BSP TTL
load( fullfile(bsp_dir, 'bsp_TTL.mat') );
bsp_TTL_ts_msec = 1e-6.*bsp_TTL_ts_ns';
bsp_TTL_ts_msec=unique(bsp_TTL_ts_msec);%for bug that I had the same TTL twice!

bsp_TTL_intervals = diff(bsp_TTL_ts_msec);
bsp_TTL_intervals_inc = diff(bsp_TTL_intervals);

%% NLG TTL
nlg_TTL_file_name = fullfile(aud_dir, 'EVENTS__Digital in.nev');
FieldSelection = [1 0 0 0 0];
ExtractHeader = 0;
ExtractMode = 1;
ModeArray = [];
nlg_TTL_ts_usec = Nlx2MatEV( nlg_TTL_file_name ,FieldSelection,ExtractHeader,ExtractMode,ModeArray);
nlg_TTL_ts_msec = nlg_TTL_ts_usec*1e-3;
nlg_TTL_intervals = diff(nlg_TTL_ts_msec);
nlg_TTL_intervals_inc = diff(nlg_TTL_intervals);

%% sync (find matching TTL intervals)
thr = 2;
x = diff(bsp_TTL_ts_msec);
y = diff(nlg_TTL_ts_msec);
[dist,ix,iy] = dtw(x,y);
pairs = [x(ix);y(iy)];
rsdl = diff(pairs);
IX = find( abs(rsdl) < thr  );
% 
% % TEST what if there are couples with more than one match:
% u=unique(ix(~isnan(ix)));
% n=histc(ix,u);
% d = u(n > 1);
% out_x_ind =find(ismember(ix,d));
% out_x_ind=union(out_x_ind,out_x_ind+1);
% 
% u=unique(iy(~isnan(iy)));
% n=histc(iy,u);
% d = u(n > 1);
% out_y_ind =find(ismember(iy,d));
% out_y_ind=union(out_y_ind,out_y_ind+1);
% 
% IX(union(out_x_ind, out_y_ind))=[];

pairs_good = pairs(:,IX);
pairs_good_bsp_ts = mean(bsp_TTL_ts_msec([ix(IX), ix(IX)+1]) .* 1e6,2);

mathing_TTL_bsp_ts = bsp_TTL_ts_msec(union(ix(IX), ix(IX)+1))*1e6;
mathing_TTL_nlg_ts = nlg_TTL_ts_msec(union(iy(IX), iy(IX)+1))*1e3;

if length(unique(pairs(1,IX))) ~= length(IX) || ...
   length(unique(pairs(2,IX))) ~= length(IX)
   error()
end

%% correct for sync jump
if ~isempty(sync_jump_ts)
    ti = [-inf sync_jump_ts inf];
    ti = [ti(1:end-1); ti(2:end)]';
    [~, IX_per_ti] = get_data_in_ti(mathing_TTL_bsp_ts, ti);

    new_ttl_bsp = [];
    new_ttl_nlg = [];
    for ii_jump = 1:length(sync_jump_ts)
        TTL_IX_before = IX_per_ti{ii_jump};
        TTL_IX_after  = IX_per_ti{ii_jump+1};
        before_jump_bsp_ts = sync_jump_ts(ii_jump) - 1;
        after_jump_bsp_ts  = sync_jump_ts(ii_jump) + 1;
        before_jump_nlg_ts = interp1(mathing_TTL_bsp_ts(TTL_IX_before), mathing_TTL_nlg_ts(TTL_IX_before), before_jump_bsp_ts, 'linear','extrap');
        after_jump_nlg_ts  = interp1(mathing_TTL_bsp_ts(TTL_IX_after ), mathing_TTL_nlg_ts(TTL_IX_after ), after_jump_bsp_ts,  'linear','extrap');
        new_ttl_bsp = [new_ttl_bsp before_jump_bsp_ts after_jump_bsp_ts];
        new_ttl_nlg = [new_ttl_nlg before_jump_nlg_ts after_jump_nlg_ts];
    end

    sync_ts.bsp = sort([mathing_TTL_bsp_ts new_ttl_bsp]);
    sync_ts.nlg = sort([mathing_TTL_nlg_ts new_ttl_nlg]);
else
    sync_ts.bsp = mathing_TTL_bsp_ts;
    sync_ts.nlg = mathing_TTL_nlg_ts;
end

%% create sync figure
ax_h = [];
figure('Units','normalized','Position',[0 0 1 1]);
pnl = panel();
pnl.pack('h',2);
pnl(1).pack('v',2);
pnl(2).pack('v',3);
pnl.margin = 30;
h=pnl.title(bsp_dir);h.Position = [0.5 1.06]; h.FontSize=16;

% plot TTL intervals
pnl(1,1).select(); hold on;
x = (bsp_TTL_ts_msec(2:end) - bsp_TTL_ts_msec(1)) .* 1e-3/60;
y = bsp_TTL_intervals*1e-3;
plot(x,y,'.b')
text(0.8,0.99, sprintf('n=%d',length(x)), 'Units','normalized','Color','b');
x = (nlg_TTL_ts_msec(2:end) - nlg_TTL_ts_msec(1)) .* 1e-3/60;
y = nlg_TTL_intervals*1e-3;
plot(x,y,'or')
text(0.8,0.95, sprintf('n=%d',length(x)), 'Units','normalized','Color','r');
xlabel('Time (minutes)')
ylabel('inter-TTL-interval (sec)')
legend({'bsp';'nlg'})
title('TTL intervals timings (relative to 1st TTL)')

pnl(1,2).select(); hold on;
x = (bsp_TTL_ts_msec(3:end) - bsp_TTL_ts_msec(1)) .* 1e-3/60;
y = bsp_TTL_intervals_inc;
plot(x,y,'.b')
text(0.8,0.99, sprintf('n=%d',length(x)), 'Units','normalized','Color','b');
x = (nlg_TTL_ts_msec(3:end) - nlg_TTL_ts_msec(1)) .* 1e-3/60;
y = nlg_TTL_intervals_inc;
plot(x,y,'or')
text(0.8,0.95, sprintf('n=%d',length(x)), 'Units','normalized','Color','r');
xlabel('Time (minutes)')
ylabel('inter-TTL-interval increament (msec)')
legend({'bsp';'nlg'})
title('TTL interval increament timings (relative to 1st TTL)')

% validate resonable gain
pnl(2,1).select(); hold on;
plot(pairs_good_bsp_ts, pairs_good(1,:)./pairs_good(2,:), '.-')
plot(repmat(sync_jump_ts,2,1), repmat(get(gca,'ylim'),length(sync_jump_ts),1)', 'm-')
rescale_plot_data('x', [1e-9/60 0]);
xlabel('bsp time (minutes)')
ylabel('clock gain')
title('only matching intervals')

pnl(2,2).select(); hold on;
plot(edges2centers(mathing_TTL_bsp_ts), 1e-3.*diff(mathing_TTL_bsp_ts) ./ diff(mathing_TTL_nlg_ts),'.-')
plot(repmat(sync_jump_ts,2,1), repmat(get(gca,'ylim'),length(sync_jump_ts),1)', 'm-')
rescale_plot_data('x', [1e-9/60 0]);
xlabel('bsp time (minutes)')
ylabel('clock gain')
title('all TTLs from matching intervals')

pnl(2,3).select(); hold on;
plot(edges2centers(sync_ts.bsp), 1e-3.*diff(sync_ts.bsp) ./ diff(sync_ts.nlg),'.-')
plot(repmat(sync_jump_ts,2,1), repmat(get(gca,'ylim'),length(sync_jump_ts),1)', 'm-')
rescale_plot_data('x', [1e-9/60 0]);
xlabel('bsp time (minutes)')
ylabel('clock gain')
title('after adding dummy TTLs in sync jumps')

linkaxes(pnl(2).de.axis,'x')

saveas(gcf, fullfile(out_dir, 'sync_bsp_nlg__TTL_timinigs'), 'jpeg')
saveas(gcf, fullfile(out_dir, 'sync_bsp_nlg__TTL_timinigs'), 'fig')

%% save interpolation values for later use
save( fullfile(out_dir, 'bsp2nlg_sync_ts') , 'sync_ts');

%% convert nlg timestamps:

 ts_synced = interp1(sync_ts.nlg, sync_ts.bsp, ts_to_convert, 'linear','extrap');

% %% convert bsp data timestamps
% files_to_convert = dir(fullfile(bsp_dir, '*pos*.mat'));
% for ii_tag = 1:length(bsp_data)
%     clear bsp_pos
%     %load( fullfile(bsp_dir, files_to_convert(ii_file).name) );
%     bsp_pos.ts_nlg_usec = interp1(sync_ts.nlg, sync_ts.nlg, bsp_data(ii_tag).ts_ns, 'linear','extrap');
%     bsp_data(ii_tag).ts_nlg_usec=bsp_pos.ts_nlg_usec;
%     %save( fullfile(bsp_dir, files_to_convert(ii_file).name), 'bsp_pos');
% end
