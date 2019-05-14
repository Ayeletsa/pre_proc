function [vsp_ts_msec_fitted_to_nlg,pp,chosen_slope] = sync_nlg_to_vsp_with_inital_fit(TTL_ts_1, TTL_ts_2, N_smoothing, max_diff, vsp_time_vec_msecs, out_dir)
% TTL ts in msec !!!
% TTL1 is vsp, TTL2 is bsp
% TTL 1 and 2 can have different number of TTLs
% p is the linear polyfit that converts time_1 to time_2
% p_inverse is the linear polyfit that converts time_2 to time_1

%%

if isempty(N_smoothing)
    N_smoothing = 5;
end


%% create time-series vector for TTL2

TTL_ts_2=round(TTL_ts_2);
TTL_vec_2 = zeros(1, TTL_ts_2(end)-TTL_ts_2(1)+1);
TTL_vec_2( ceil(TTL_ts_2-TTL_ts_2(1)) + 1 ) = 1;
zero_padd_duration = 60*1e3;
zero_padd = zeros(1,zero_padd_duration);
TTL_vec_2 = [zero_padd TTL_vec_2 zero_padd];
if N_smoothing > 1
    TTL_vec_2 = filtfilt(hamming(N_smoothing),1,TTL_vec_2);
end



%% find initial slope and cross correlate (find shift)

slopes=(0.99:0.001:1.01);
max_corr=zeros(length(slopes),1);
IX_max_corr=zeros(length(slopes),1);
TTL_vec_1_array=cell(length(slopes),1);
corr_array=cell(length(slopes),1);
lags_array=cell(length(slopes),1);

for ii_slope=1:length(slopes)
    TTL_ts_1_modified = round (TTL_ts_1 * slopes(ii_slope));
    TTL_vec_1 = zeros(1, TTL_ts_1_modified(end)-TTL_ts_1_modified(1)+1);
    TTL_vec_1( ceil(TTL_ts_1_modified-TTL_ts_1_modified(1)) + 1 ) = 1;
    TTL_vec_1 = [zero_padd TTL_vec_1 zero_padd];
    if N_smoothing > 1
        TTL_vec_1_array{ii_slope} = filtfilt(hamming(N_smoothing),1,TTL_vec_1);
    end
    
    [corr_array{ii_slope} lags_array{ii_slope}] = xcorr(TTL_vec_1_array{ii_slope}, TTL_vec_2);
    [max_corr(ii_slope), IX_max_corr(ii_slope)] = max(corr_array{ii_slope});
    
end

[~,IX] = max(max_corr);
chosen_slope = slopes(IX);
TTL_ts_1_sloped = round (TTL_ts_1 * chosen_slope);
c=corr_array{IX};
lags=lags_array{IX};
IX_lags=IX_max_corr(IX);
shift_msec = -lags(IX_lags) + TTL_ts_2(1) - TTL_ts_1_sloped(1);
TTL_ts_1_shifted = TTL_ts_1_sloped + shift_msec;

[PKS LOCS] = findpeaks(c,'MINPEAKHEIGHT',max(c)/100);
[PKS_sorted PKS_sorted_IX] = sort(PKS, 'descend');
LOCS_sorted = LOCS(PKS_sorted_IX);
sidelobe_ratio = PKS_sorted(1) / PKS_sorted(2);
sidelobe_dB = 20*log10(sidelobe_ratio);
sidelobe_delay = abs(LOCS_sorted(1) - LOCS_sorted(2));



%% find best TTL pairs

temp_diff = [];
IX_closest_TTL_2 = [];
number_of_pairs = min( length(TTL_ts_1_sloped), length(TTL_ts_2) );
for ii_TTL = 1:length(TTL_ts_1_shifted)
    [temp_diff(ii_TTL), IX_closest_TTL_2(ii_TTL)] = min( abs(TTL_ts_1_shifted(ii_TTL) - TTL_ts_2) );
end
[temp_diff_sorted,  IX_TTL_1_temp_diff_sorted] = sort(temp_diff,'ascend');
IX_matching_TTL_1 = IX_TTL_1_temp_diff_sorted(1:number_of_pairs);
IX_matching_TTL_2 = IX_closest_TTL_2(IX_matching_TTL_1);
IX_matching_TTL_1 = sort(IX_matching_TTL_1);
IX_matching_TTL_2 = sort(IX_matching_TTL_2);


%% remove pairs with temporal diff > max_diff

X = TTL_ts_1_shifted(IX_matching_TTL_1);
Y = TTL_ts_2(IX_matching_TTL_2);
bad_pairs_IX = find( abs(X-Y) > max_diff);
IX_matching_TTL_1(bad_pairs_IX) = [];
IX_matching_TTL_2(bad_pairs_IX) = [];


%% 1st fit (before removing bad points)

X = TTL_ts_1_sloped(IX_matching_TTL_1);
Y = TTL_ts_2(IX_matching_TTL_2);
vsp_p = polyfit(X,Y,1);
vsp_p_inverse = polyfit(Y,X,1);


%% remove pairs with temporal diff > max_diff

X = TTL_ts_1_sloped(IX_matching_TTL_1);
Y = TTL_ts_2(IX_matching_TTL_2);
bad_pairs_IX = find( abs(polyval(vsp_p,X)-Y) > max_diff);
IX_matching_TTL_1(bad_pairs_IX) = [];
IX_matching_TTL_2(bad_pairs_IX) = [];


%% 2nd fit (after removing bad points)

X = TTL_ts_1_sloped(IX_matching_TTL_1);
Y = TTL_ts_2(IX_matching_TTL_2);
vsp_p = polyfit(X,Y,1);
vsp_p_inverse = polyfit(Y,X,1);

%% use spline to create pairwise fit

pp = spline(X,Y);

pp_odd = spline(X(1:2:end),Y(1:2:end));
xx_even = ppval(pp_odd,X(2:2:end));

pp_even = spline(X(2:2:end),Y(2:2:end));
xx_odd = ppval(pp_even,X(1:2:end));

%% plot figure

figure('units','normalized','outerposition',[0 0 1 1],'color','w')

subplot(4,4,1:4)
vsp_tone_times_intervals = diff(TTL_ts_1);
nlx_tone_intervals = diff(TTL_ts_2);
plot(vsp_tone_times_intervals)
hold on
plot(nlx_tone_intervals)
legend('Audiologger Tone Intervals', 'nlx Tone Intervals')
n_beeps = max(length(vsp_tone_times_intervals),length(nlx_tone_intervals));
xlim([0 n_beeps+1])
xlabel('# beep')
ylabel('time (ms)')
title('Inter beep interval')
box off

subplot(4,4,5)
hold on
plot (max_corr,'marker','*')
s_h=scatter(IX,max_corr(IX),[],'r','marker','*');
set(gca,'xtick',1:3:length(slopes),'xticklabel',slopes(1:3:end));
legend(s_h,['Chosen slope = ' num2str(chosen_slope)]);
xlabel('slope');
ylabel('max corr');
xlim([0 length(slopes)+1])
title('initial slope')

subplot(4,4,6)
hold on
plot(lags,c)
plot(lags(LOCS),c(LOCS),'or')
text(lags(LOCS_sorted(2))+200,c(LOCS_sorted(2)), ...
    {['1^{st} sidelobe: ' num2str(sidelobe_dB) 'dB'];...
    [num2str(sidelobe_delay) 'ms lag']},...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
xlabel('lags')
ylabel('xcorr')

subplot(4,4,7)
hold on
plot(X,Y,'o')
plot(X,polyval(vsp_p,X),'r')
legend({'data points';'polyval'})
title(['1^{st} fit (before removing bad points), n=' num2str(length(X))])

subplot(4,4,8)
hold on
plot(X,Y,'o')
plot(X,polyval(vsp_p,X),'r')
legend({'data points';'polyval'})
title(['2^{nd} fit (after removing bad points), n=' num2str(length(X))])

subplot(4,4,9:12)
plot(1:length(Y),polyval(vsp_p,X)-Y,'marker','*')
box off
xlim([0 length(Y)+1])
% xlabel('# beep')
ylabel('Time (ms)')
legend('fitted vsp - nlg')
title('Time differnce between beeps after 2^{nd} fit')

subplot(4,4,13:16)
hold on
plot([1:2:length(X)],xx_odd - Y(1:2:end),'color',[.8 0 .8])
plot([2:2:length(X)],xx_even - Y(2:2:end),'color',[0 .8 .8])
box off
xlim([0 length(Y)+1])
xlabel('# beep')
ylabel('Time (ms)')
legend('Odd couples','Even couples')
title('Time differnce after pair-wise fit on half of the data')

%% 

saveas(gcf, fullfile(out_dir, 'sync_vsp2nlg'), 'jpeg')
saveas(gcf, fullfile(out_dir, 'sync_vsp2nlg'), 'fig')
close(gcf)
save(fullfile(out_dir, 'time_conv_pp_msec_vsp2nlg'), 'pp')
save(fullfile(out_dir, 'slope_for_vsp_ts'), 'chosen_slope')

vsp_ts_msec_fitted_to_nlg = ppval(pp,vsp_time_vec_msecs * chosen_slope);


end


