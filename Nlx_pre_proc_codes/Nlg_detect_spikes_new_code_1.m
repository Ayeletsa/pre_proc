function p =Nlg_detect_spikes_new_code (p)
main_dir = p.path_day_dir;
active_TTs = p.use_tetrodes;
active_channels = p.active_channels;
%% New Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function Nlg_detect_spikes_CSC_2_1(expname,forcerecalc)
%
% based on Didi Omer(November 2015) modified by Ayelet Sarel (Oct 2017)
%
%
date = char (regexp (main_dir,'\d{8}','match'));
main_dir_out = [p.path_dataout,'\',p.year_bat_path,'\',date];
csc_dir_spikes = fullfile(main_dir_out,'spikes');
base_output_dir = fullfile(main_dir_out,'detected_spikes');

if exist (base_output_dir,'dir')
    return
else
    mkdir(base_output_dir)
end

bat_id = num2str(p.bat);

%%
pool = gcp('nocreate'); % If no pool, do not create new one.
if isempty(pool)
    poolsize = 0;
else
    poolsize = pool.NumWorkers;
end
if poolsize ~=6  % checking to see if my pool is already open
    parpool(6)
else
    disp('matlab workers already open')
end

% Parameters
% TODO: call them from p.in and difeine them with right names and params
%-----------------------------------------------------------------
include_negative_threshold=p.include_negative_threshold; %find also negative events
do_coincidence_detection=p.do_coincidence_detection;
compare_to_library=p.compare_to_library;

threhold_factor=p.threhold_factor; %for voltage threshold (threhold=median*factor)
min_sep_events=p.min_sep_events; %min separation between events - ~=500us
r_threshold=p.r_threshold;
coincidence_window=p.coincidence_window; % Duration of coincidence-detection window (in microsec) ACROSS tetrodes
thresh_for_coinc=p.thresh_for_coinc; %remove if X tetrodes has event in the window

create_movie=p.create_movie;

% include_negative_threshold=1; %find also negative events
% do_coincidence_detection=1;
% compare_to_library=1;
%
% threhold_factor=4.6; %for voltage threshold (threhold=median*factor)
% min_sep_events=15; %min separation between events - ~=500us
% r_threshold = 0.9;
% coincidence_window=500; % Duration of coincidence-detection window (in microsec) ACROSS tetrodes
% thresh_for_coinc=4; %remove if X tetrodes has event in the window
%
% create_movie=0;

%%%%% remove the following parameters in your codes if you have them in a
%%%%% different format!!
% active_TTs=1:4; %how many tetrodes are connected for coincidence detection
% active_channels=[0,1,0,1;0,1,1,1;1,1,1,0; 0,0,1,1];

tic1 = tic;
disp('start: Nlg_detect_spikes_CSC_3_6()');

library_file_name=p.library_file_name;

load(library_file_name); % Load the library of acceptable spike shapes;
library_of_acceptable_spike_shapes=new_lib;

% forcerecalc=1;
% folder_name = 'D:\Matlab\Spike_detection_proj\Data_for_library\yr2014_bat4500_Neo_Nlg\spike_detection_20140921';
% if exist(folder_name,'dir') && ~forcerecalc
%     disp(['Detect-spikes was already done']);
%     return;
% end


thresholds = zeros(4,4);

for ii_tetrode=active_TTs
    
    tic2 = tic;
    filename = cell(1,4);
    Spiketimeidx = cell(1,4);
    csc = cell(4,1);
    %path = [csc_dir_spikes,'spikes__TT'];
    %tetrodes_use_for_sorting = param.tetrodes.use_for_sorting;
    thresh = zeros(1,4);
    neural_data_ts=[];
    active_ch_for_tt=find(active_channels(ii_tetrode,:));
    % detect spikes
    for ch=active_ch_for_tt
        
        filename=  fullfile(csc_dir_spikes,['spikes_bat_',bat_id,'_day_',date,'_TT' num2str(ii_tetrode) '_ch' num2str(ch) '.ncs']);% originally - num2str(ch-1)
        [csc{ch}, timestamps{ch},fs] = Nlx_csc_read(filename,neural_data_ts);
        timestamps{ch}=timestamps{ch};
        if ismember(ch,active_ch_for_tt)
            csc{ch}=csc{ch};
        else
            csc{ch}=zeros(size(timestamps{ch})); %if channel is disconnected set to zero.
        end
        %         load(fullfile(csc_dir_spikes,['spikes_',bat_id,'_d',date,'_TT',num2str(ii_tetrode),'_ch',num2str(ch),'.mat']));
        %find threhols for ch:
        thresh(ch) = threhold_factor* median(abs(csc{ch}(:))); % thresh by Quiroga 2004 (Didi Jan 2017)
        %detecet positive events:
        [spikeidx_max,spikepk_max]  = peakseek(csc{ch},1,thresh(ch));
        [spikeidx_max,ii] = sort(spikeidx_max);
        spikepk_max = spikepk_max(ii);
        spikeidx_max = spikeidx_max((spikeidx_max>7 & spikeidx_max<=length(csc{ch})-24));
        spiketime_max = timestamps{ch}(spikeidx_max);
        
        
        if include_negative_threshold
            %detect negative events:
            [spikeidx_min,spikepk_min]  = peakseek(-csc{ch},1,thresh(ch));
            [spikeidx_min,ii] = sort(spikeidx_min);
            spikepk_min = spikepk_min(ii);
            spikeidx_min = spikeidx_min((spikeidx_min>7 & spikeidx_min<=length(csc{ch})-24));
            spiketime_min = timestamps{ch}(spikeidx_min);
            %unite positive and negative events:
            [spiketime,unite_ii]=sort([spiketime_max(:)',spiketime_min(:)']);
            spikeidx=[spikeidx_max(:)',spikeidx_min(:)'];
            spikeidx=spikeidx(unite_ii);
            spiketime_all_ch{ch,1}=spiketime;
            spikeidx_all_ch{ch,1}=spikeidx;
        else
            [spiketime,unite_ii]=sort([spiketime_max(:)']);
            spikeidx=[spikeidx_max(:)'];
            spikeidx=spikeidx(unite_ii);
            spiketime_all_ch{ch,1}=spiketime;
            spikeidx_all_ch{ch,1}=spikeidx;
        end
    end
    
    % TODO: save treshold in P after running on all TTs
    threshold(:,ii_tetrode)= thresh(:);
    %     SamplingFreq=Spikes.params.CSC_SamplingFreq;
    SamplingFreq=fs;
    disp(['1) detecting spikes on TT',num2str(ii_tetrode),' elepsed:',num2str(toc(tic2)/60),' min']);
    
    %% unite events from all channels
    tic3 = tic;
    startidx = min([spikeidx_all_ch{1}' ;spikeidx_all_ch{2}';spikeidx_all_ch{3}';spikeidx_all_ch{4}']);
    endidx  =  max([spikeidx_all_ch{1}' ;spikeidx_all_ch{2}';spikeidx_all_ch{3}';spikeidx_all_ch{4}']);
    
    
    len = endidx-startidx  +1;
    M = zeros(4,len);
    for i=active_ch_for_tt
        M(i,spikeidx_all_ch{i})=csc{i,1}(spikeidx_all_ch{i});
    end
    M = max(abs(M));
    
    spikeidx = peakseek(M,min_sep_events,min(thresh(:)));
    spiketime = timestamps{ch}(spikeidx);
    clear M;
    clear xx;
    clear idx;
    disp(['2) remove redundent spike events on TT',num2str(ii_tetrode),' elepsed:',num2str(toc(tic3)/60),' min']);
    
    %% calculate spikewave
    tic4=tic;
    spikewave = zeros(32,size(spiketime,2),4);
    mat= repmat([-7:24]',1,length(spikeidx));
    spike_idx = repmat(spikeidx,[32,1]);
    spike_idx = spike_idx+mat;
    
    for ch=active_ch_for_tt
        spikewave(:,:,ch) = csc{ch}(spike_idx);
    end
    
    disp(['3) calculating waveforms TT',num2str(ii_tetrode),' (',num2str(size(spiketime,2)),' spikes)',' elepsed:',num2str(toc(tic4)/60),' min']);
    clear csc;
    
    if create_movie
        video_file_name=[base_output_dir,'/threshold_mov_TT',num2str(ii_tetrode)];
        create_thresh_movie(video_file_name,timestamps,spiketime,spikewave,thresh)
    end
    
    %% compare to acceptable spike shapes
    if compare_to_library==1
        tic5 = tic;
        
        [~,iimax ] = max(squeeze(max(abs(spikewave))),[],2);
        waves = zeros(size(spikewave,1),size(spikewave,2));
        for i=1:size(waves,2)
            waves(:,i) = spikewave(:,i,iimax(i));
        end
        
        c = zeros(size(library_of_acceptable_spike_shapes,2),size(waves,2));
        
        parfor ii=1:size(library_of_acceptable_spike_shapes,2)
            mat = repmat(library_of_acceptable_spike_shapes(:,ii),1,size(waves,2));
            c1 = corrcoeff(waves(2:end-1,:),mat(2:end-1,:));
            c2 = corrcoeff(waves(1:end-2,:),mat(2:end-1,:));
            c3 = corrcoeff(waves(3:end,:),mat(2:end-1,:));
            c(ii,:) = max([c1;c2;c3]);
        end
        
        
        
        c = max(c);
        acceptedidx = find(c>=r_threshold);
        spike_idx_TT{ii_tetrode}=spikeidx(acceptedidx);
        Timestamps_accepted_spikes_TT{ii_tetrode}=spiketime(acceptedidx);
        spikes_TT{ii_tetrode}=spikewave(:,acceptedidx,:);
        disp(['4) Testing spikewaves against spike-shape lib. TT',num2str(ii_tetrode),' elepsed:',num2str(toc(tic5)/60),' min' ]);
        
        
        
        figure
        hist(c,1000)
        hold on
        plot([0.9 0.9], [0 max(hist(c,1000))],'r')
        title('max r to library')
        fig_name=[base_output_dir,'library_r_TT',num2str(ii_tetrode),'.jpg'];
        saveas(gcf,fig_name)
        close(gcf)
        clear c;
        
    else
        
        Timestamps_accepted_spikes_TT{ii_tetrode}=spiketime(:)';
        spikes_TT{ii_tetrode}=spikewave(:,:,:);
        
        
        
        
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Coincidence-Detection across Tetrodes to eliminate artefacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic8 = tic;
if length(active_TTs)>1
    idx_coincidence_vec = cell(length(active_TTs),1);
    
    for ii_spikes = 1:length(Timestamps_accepted_spikes_TT{1}) % Loop over the spikes of the FIRST tetrode
        temp_stack = cell(length(active_TTs),1);
        t_spike = Timestamps_accepted_spikes_TT{1}(ii_spikes);
        temp_stack{1,1} = ii_spikes;
        for TT=2:length(active_TTs)
            temp_stack{TT,1} = find( abs( Timestamps_accepted_spikes_TT{TT} - t_spike ) <= coincidence_window ); % THE COINCIDENCE DETECTION
        end
        ii = ~cellfun(@isempty,temp_stack);
        if sum(ii)>=length(active_TTs)
            for jj=1:length(active_TTs)
                idx_coincidence_vec{jj,1} = cat(2,idx_coincidence_vec{jj,1}, temp_stack{jj,1});
            end
        end
    end
end

% %
for i=1:length(active_TTs)
    idx_coincidence_vec{i} = unique(idx_coincidence_vec{i});
    if ~isempty(idx_coincidence_vec{i})
        Timestamps_accepted_spikes_TT{i}(idx_coincidence_vec{i})=[];
        spikes_TT{i}(:,idx_coincidence_vec{i},:)=[];
    end
end
disp(['5) cross TT coincident spike removal...',' elepsed:',num2str(toc(tic8)/60),' min']);



%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the data in Neuralynx NTT files for three cases:
%%%%%%%%%%%%%%%%%%%%%%%%%

tic9 = tic;
for ii_Tetrode=1: length(active_TTs)
    
    %base_output_dir						= [param.path.analysis, 'preprocessing_spikes4',filesep];
    if ~exist(base_output_dir,'dir'), mkdir(base_output_dir); end
    filename_out = [base_output_dir, '\spikes__TT', num2str(active_TTs(ii_Tetrode))];
    
    % load generic header to add to NTT files
    load('ntt_generic_header.mat');
    
    %     l = ['-ThreshVal ',num2str(threshold(1,ii_Tetrode)),' ',...
    %                     num2str(threshold(2,ii_Tetrode)),' ',...
    %                     num2str(threshold(3,ii_Tetrode)),' ',...
    %                     num2str(threshold(4,ii_Tetrode))];
    %     Header{end+1} = l;
    %     Header = MakeNlxHeader(Header);
    %Header = Header';
    
    Timestamps_accepted_spikes=Timestamps_accepted_spikes_TT{ii_Tetrode};
    spikes=spikes_TT{ii_Tetrode};
    spikes = permute(spikes,[1,3,2]);
    
    %only export timestamps and data points - Full session
    FieldSelection = [1 0 0 0 1 1];
    
    if ispc
        Mat2NlxSpike([filename_out,'.NTT'],0,1,[],FieldSelection,Timestamps_accepted_spikes(:)', spikes,Header);
    elseif isunix
        if exist([filename_out,'.NTT'],'file'), delete([filename_out,'.NTT']),end
        Mat2NlxTT([filename_out,'.NTT'],0,1,1,length(Timestamps_accepted_spikes(:)),FieldSelection,Timestamps_accepted_spikes(:)', spikes,Header);
    end
end
disp(['6) Saving spikes to disk...',' elepsed:',num2str(toc(tic8)/60),' min']);
disp(['8) Total detection on all TTs elepsed:',num2str(toc(tic1)/60),' min']);

p.voltage_treshold_for_detection = threshold;






