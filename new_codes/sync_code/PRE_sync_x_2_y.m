function sync_ts=PRE_sync_x_2_y(TTL_x,TTL_y,x_name,y_name,sync_name,sync_dir)


X_TTL_intervals = diff(TTL_x);
X_TTL_intervals_inc = diff(X_TTL_intervals);


Y_TTL_intervals = diff(TTL_y);
Y_TTL_intervals_inc = diff(Y_TTL_intervals);

%% sync (find matching TTL intervals)
thr = 2;
x = diff(TTL_x);
y = diff(TTL_y);
[dist,ix,iy] = dtw(x,y);
pairs = [x(ix);y(iy)];
rsdl = diff(pairs);
IX = find( abs(rsdl) < thr  );

pairs_good = pairs(:,IX);
pairs_good_x_ts = mean(TTL_x([ix(IX), ix(IX)+1]) .* 1e6,2);

mathing_TTL_X_ts = TTL_x(union(ix(IX), ix(IX)+1))*1e6;
mathing_TTL_Y_ts = TTL_y(union(iy(IX), iy(IX)+1))*1e3;

if length(unique(pairs(1,IX))) ~= length(IX) || ...
   length(unique(pairs(2,IX))) ~= length(IX)
   error()
end

%% correct for sync jump
% TODO add this option if there will be bugs!
sync_jump_ts=[];

if ~isempty(sync_jump_ts)
    ti = [-inf sync_jump_ts inf];
    ti = [ti(1:end-1); ti(2:end)]';
    [~, IX_per_ti] = get_data_in_ti(mathing_TTL_X_ts, ti);

    new_ttl_X = [];
    new_ttl_Y = [];
    for ii_jump = 1:length(sync_jump_ts)
        TTL_IX_before = IX_per_ti{ii_jump};
        TTL_IX_after  = IX_per_ti{ii_jump+1};
        before_jump_X_ts = sync_jump_ts(ii_jump) - 1;
        after_jump_X_ts  = sync_jump_ts(ii_jump) + 1;
        before_jump_Y_ts = interp1(mathing_TTL_X_ts(TTL_IX_before), mathing_TTL_Y_ts(TTL_IX_before), before_jump_X_ts, 'linear','extrap');
        after_jump_Y_ts  = interp1(mathing_TTL_X_ts(TTL_IX_after ), mathing_TTL_Y_ts(TTL_IX_after ), after_jump_X_ts,  'linear','extrap');
        new_ttl_X = [new_ttl_X before_jump_X_ts after_jump_X_ts];
        new_ttl_Y = [new_ttl_Y before_jump_Y_ts after_jump_Y_ts];
    end

    sync_ts.X = sort([mathing_TTL_X_ts new_ttl_X]);
    sync_ts.Y = sort([mathing_TTL_Y_ts new_ttl_Y]);
else
    sync_ts.X = mathing_TTL_X_ts;
    sync_ts.Y = mathing_TTL_Y_ts;
end

%% create sync figure
ax_h = [];
figure('Units','normalized','Position',[0 0 1 1]);
pnl = panel();
pnl.pack('h',2);
pnl(1).pack('v',2);
pnl(2).pack('v',3);
pnl.margin = 30;
h=pnl.title(x_name);h.Position = [0.5 1.06]; h.FontSize=16;

% plot TTL intervals
pnl(1,1).select(); hold on;
x = (TTL_x(2:end) - TTL_x(1)) .* 1e-3/60;
y = X_TTL_intervals*1e-3;
plot(x,y,'.b')
text(0.8,0.99, sprintf('n=%d',length(x)), 'Units','normalized','Color','b');
x = (TTL_y(2:end) - TTL_y(1)) .* 1e-3/60;
y = Y_TTL_intervals*1e-3;
plot(x,y,'or')
text(0.8,0.95, sprintf('n=%d',length(x)), 'Units','normalized','Color','r');
xlabel('Time (minutes)')
ylabel('inter-TTL-interval (sec)')
legend({x_name;y_name})
title('TTL intervals timings (relative to 1st TTL)')

pnl(1,2).select(); hold on;
x = (TTL_x(3:end) - TTL_x(1)) .* 1e-3/60;
y = X_TTL_intervals_inc;
plot(x,y,'.b')
text(0.8,0.99, sprintf('n=%d',length(x)), 'Units','normalized','Color','b');
x = (TTL_y(3:end) - TTL_y(1)) .* 1e-3/60;
y = Y_TTL_intervals_inc;
plot(x,y,'or')
text(0.8,0.95, sprintf('n=%d',length(x)), 'Units','normalized','Color','r');
xlabel('Time (minutes)')
ylabel('inter-TTL-interval increament (msec)')
legend({x_name;y_name})
title('TTL interval increament timings (relative to 1st TTL)')

% validate resonable gain
pnl(2,1).select(); hold on;
plot(pairs_good_x_ts, pairs_good(1,:)./pairs_good(2,:), '.-')
plot(repmat(sync_jump_ts,2,1), repmat(get(gca,'ylim'),length(sync_jump_ts),1)', 'm-')
rescale_plot_data(x_name, [1e-9/60 0]);
xlabel([x_name,' time (minutes)'])
ylabel('clock gain')
title('only matching intervals')

pnl(2,2).select(); hold on;
plot(edges2centers(mathing_TTL_X_ts), 1e-3.*diff(mathing_TTL_X_ts) ./ diff(mathing_TTL_Y_ts),'.-')
plot(repmat(sync_jump_ts,2,1), repmat(get(gca,'ylim'),length(sync_jump_ts),1)', 'm-')
rescale_plot_data(x_name, [1e-9/60 0]);
xlabel([ x_name,'time (minutes)'])
ylabel('clock gain')
title('all TTLs from matching intervals')

pnl(2,3).select(); hold on;
plot(edges2centers(sync_ts.X), 1e-3.*diff(sync_ts.X) ./ diff(sync_ts.Y),'.-')
plot(repmat(sync_jump_ts,2,1), repmat(get(gca,'ylim'),length(sync_jump_ts),1)', 'm-')
rescale_plot_data(x_name, [1e-9/60 0]);
xlabel('X time (minutes)')
ylabel('clock gain')
title('after adding dummy TTLs in sync jumps')

linkaxes(pnl(2).de.axis,'x')

saveas(gcf, fullfile(sync_dir, [sync_name,'__TTL_timinigs']), 'jpeg')
saveas(gcf, fullfile(sync_dir, [sync_name,'__TTL_timinigs']), 'fig')

%% save interpolation values for later use
save( fullfile(sync_dir, ['sync_',x_name,'2',y_name]) , 'sync_ts');

%%

