function [vsp_ts_msecs, vsp_tone_ts_msecs] = Find_Tone_Times_in_Audio (sigAudioUnfilt,main_dir)

sync_dir = fullfile(main_dir, 'sync');
out_dir = fullfile(sync_dir, 'sync_vsp2bsp');

%% input args

fsAudio = 100000;
beep_freq=4*1e3;
% max_tone_length_secs = 0.305; % harsh parameters
% min_tone_length_secs = 0.208; 
max_tone_length_secs = 0.320; % soft parameters
min_tone_length_secs = 0.280;


%% load audio signal

% sigAudioUnfilt = load(fullfile(vsp_dir,'sigAudioUnfilt'));
% sigAudioUnfilt = sigAudioUnfilt.sigAudioUnfilt;
if ~exist(out_dir,'dir')
    mkdir(out_dir)
end

%% estimate the real frequency of the recorded beeps, peaked around expected frequancy

% [tone_power,freq_pwer] = pwelch(sigAudioUnfilt,[],[],[],fsAudio);
% IX_beep_freq = find(freq_pwer > beep_freq-50 & freq_pwer < beep_freq+50);
% [~,IX] = max(tone_power(IX_beep_freq));
% estimated_beep_freq = freq_pwer(IX_beep_freq(IX));
estimated_beep_freq = beep_freq;


%% design bandpass FIR filter for tone detection

% % normalized cutoff frequencies
cutoff_low_freq = (estimated_beep_freq - 7) / (fsAudio / 2);
cutoff_high_freq = (estimated_beep_freq + 7) / (fsAudio / 2);

filt_order = 300;
% filt_delay = filt_order / 2;

% get the filter
fir_tone = fir1(filt_order, [cutoff_low_freq, cutoff_high_freq]);


%% filter the signal for tone detection

% using filtfilt
tic
sigAudio = cell(length(sigAudioUnfilt),1);
for i = 1:length(sigAudioUnfilt)
    sigAudio{i} = filtfilt(fir_tone, 1, sigAudioUnfilt{i});
end
toc


%% split signal into chunks of 50 so that hilbert doesn't run out of memory

r = 50;
m = numel(sigAudio);
sigAudio = [sigAudio;cell(ceil(m/r)*r-m,1)];
sigAudio = reshape(sigAudio,r,[]);

sigAudioUnfilt = [sigAudioUnfilt;cell(ceil(m/r)*r-m,1)];
sigAudioUnfilt = reshape(sigAudioUnfilt,r,[]);
% sigVec1=cell2mat(reshape(sigAudioUnfilt,[],1));

env = cell(1, size(sigAudio,2));
for i = 1:size(sigAudio,2)
    sigMat = cell2mat(sigAudio(:,i));
    env{i} = abs(hilbert(sigMat));
end

envelope = cell2mat(env');
sigVec=cell2mat(reshape(sigAudio,[],1));


%% find tone starts in envelope using "clustering" method
% TODO check th
% th = .5 * prctile(envelope,99.9);
th = 8e-3;

idx = find( envelope > th );
vals = envelope(idx);
idxMS = idx / fsAudio * 1000;

% Clustering
clear clusters;
clusters = struct('start', {}, 'end', {}, 'peak_abs_idx', {});
ITIs = [inf; diff(idx); inf];
ITIs_sec = ITIs / fsAudio;

% find clusters limits
% TODO limit clusters by diff from previous th crossing
% clusters_limits = find(ITIs_sec > max_gap_secs);
clusters_limits = find(ITIs > 1);
if isempty(idx)
    clusters_start = [];
    clusters_end = [];
else
    clusters_start = idx( clusters_limits(1:end-1) );
    clusters_end = idx( clusters_limits(2:end) - 1 );
    clusters( length(clusters_start) ).start = 0; % only to init array size
    clusters( length(clusters_end) ).end = 0;
end

clusters_start_cells = num2cell(clusters_start);
clusters_end_cells = num2cell(clusters_end);
[clusters(:).start] = clusters_start_cells{:};
[clusters(:).end] = clusters_end_cells{:};

% filter by cluster length
cluster_lengths_all = clusters_end - clusters_start;
clusters_idx = find(cluster_lengths_all < max_tone_length_secs * fsAudio & cluster_lengths_all > min_tone_length_secs * fsAudio);

tone_times_samples = clusters_start(clusters_idx);
env_vals = envelope(tone_times_samples);

vsp_time_vec_secs = (1:length(sigVec)) / fsAudio;
vsp_ts_msecs = vsp_time_vec_secs' * 1e3;
vsp_tone_ts_secs = tone_times_samples / fsAudio;
vsp_tone_ts_msecs = vsp_tone_ts_secs' * 1e3;



%% plot a figure of the detected beeps

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
    
%     subplot(2,3,1:2)
%     plot(vsp_time_vec_secs,sigVec)
%     hold on
%     plot(vsp_time_vec_secs,envelope)
%     plot(vsp_tone_ts_secs,env_vals,'*');
%     box off
%     xlabel ('Time (s)');
%     ylabel ('Amplitude (AU)');

tone_ITI = diff(vsp_tone_ts_secs);
n_missed_tones = sum(floor(tone_ITI / 30) - 1);
n_identified_tones = length(vsp_tone_ts_secs);

subplot(2,3,1:3)
hold on
plot(tone_ITI,'color',[.6 0 .8],'marker','.');
x_vec = 1:length(tone_ITI);
high_ITI_ind = (floor(tone_ITI / 30) - 1) > 0;
scatter(x_vec(high_ITI_ind),tone_ITI(high_ITI_ind),[],[1 .5 0],'marker','.');
box off
xlabel('# Tone')
ylabel ('ITI (s)')

subplot(2,3,4)
histogram(tone_ITI,'FaceColor',[.6 0 .8])
box off
title ('ITI''s variability')
xlabel('ITI (s)')
ylabel('# Counts')

subplot(2,3,5)
label1 = [num2str(n_identified_tones) ' Identified'];
label2 = [num2str(n_missed_tones) ' Missed'];
pie([length(vsp_tone_ts_secs) n_missed_tones],{label1,label2})
colormap([.6 0 .8;1 .5 0]);  
title ('Proportions')

subplot(2,3,6)
time_before_after_ms=100;
samples_before_after=fsAudio*time_before_after_ms/1e3;
n_beeps=length(vsp_tone_ts_msecs);
hold on

for ii_beep = 1:n_beeps
    IX_beep=(tone_times_samples(ii_beep)-samples_before_after):(tone_times_samples(ii_beep)+samples_before_after);
    h = plot(vsp_ts_msecs(IX_beep) - vsp_tone_ts_msecs(ii_beep),envelope(IX_beep));
    c = get(h,'color');
    plot(0,env_vals(ii_beep),'*','color',c);
end

axis tight
title('Tone identification')
xlabel ('Relative time to identification (ms)');
ylabel ('Amplitude (AU)');

date = char(regexp(out_dir,'\d{8}','match'));
suptitle(['Identifying tones in audio recordings - ' date])

saveas(gcf, fullfile(out_dir, 'identify_tones'), 'jpeg')
saveas(gcf, fullfile(out_dir, 'identify_tones'), 'fig')
close(gcf)

end