function p =Nlg_detect_spikes_library_first_findpeaks_1(p);
main_dir = p.path_day_dir;
active_TTs = p.use_tetrodes;
active_channels = p.active_channels;

%% New Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function Nlg_detect_spikes_CSC_2_1(expname,forcerecalc)
%
% based on Didi Omer(November 2015) modified by Ayelet Sarel (Oct 2017)
% modified by Ayelet Sarel (Oct 2017) main modifications:
% 1. Changed the function peakseek to findpeaks--> slower but detect more
%    peaks when there are few peaks in a sequnce
% 2. Changed the order- threshold-->library-->window of separation between
%     spikes-->coincidence
% 3. Added the possibility to detect negative events (for inverted spikes)
% 
%

date = char (regexp (main_dir,'\d{8}','match'));
main_dir_out = [p.path_dataout,'\',p.year_bat_path,'\',date];
csc_dir_spikes = fullfile(main_dir_out,'spikes');
base_output_dir = fullfile(main_dir_out,'\detected_spikes\');

if exist (base_output_dir,'dir')
    return
else
    mkdir(base_output_dir)
end

bat_id = num2str(p.bat);
%%

pp = gcp('nocreate'); % If no pool, do not create new one.
if isempty(pp)
    poolsize = 0;
else
    poolsize = pp.NumWorkers;
end
if poolsize ~=6  % checking to see if my pool is already open
    parpool(6)
else
    disp('matlab workers already open')
end

tic1 = tic;
disp('start: Nlg_detect_spikes_CSC_3_6()');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% remove the following parameters in your codes if you have them in a
%%%%% different format!!

% Parameters
%-----------------------------------------------------------------
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
use_for_sorting=p.use_for_sorting;
create_movie=p.create_movie;

library_file_name=p.library_file_name;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load(library_file_name); % Load the library of acceptable spike shapes;
library_of_acceptable_spike_shapes=new_lib;
thresholds = zeros(4,4);

for ii_tetrode=use_for_sorting
    
    tic2 = tic;
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
        [csc{ch}, timestamps{ch}] = Nlx_csc_read(filename,neural_data_ts);
        timestamps{ch}=timestamps{ch};
        if ismember(ch,active_ch_for_tt)
            csc{ch}=csc{ch};
        else
            csc{ch}=zeros(size(timestamps{ch})); %if channel is disconnected set to zero.
        end
        %in the following file we have some parameters for later (created in the
        %filtering code
        %find threhols for ch:
        thresh(ch) = threhold_factor* median(abs(csc{ch}(:))); % thresh by Quiroga 2004 (Didi Jan 2017)
        %detecet positive events:
        [spikepk_max,spikeidx_max]  = findpeaks(csc{ch},'MinPeakHeight',thresh(ch));
        %remove spikes from beginning or end of data:
        spikeidx_max = spikeidx_max((spikeidx_max>7 & spikeidx_max<=length(csc{ch})-24));
        spiketime_max = timestamps{ch}(spikeidx_max);
        
       
        if include_negative_threshold
            %detect negative events:
            [spikepk_min,spikeidx_min]  = findpeaks(-csc{ch},'MinPeakHeight',thresh(ch));         
             %remove spikes from beginning or end of data:
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
    
    
    thresholds(:,ii_tetrode)= thresh(:);
    disp(['1) detecting spikes on TT',num2str(ii_tetrode),' elepsed:',num2str(toc(tic2)/60),' min']);
    
    %% unite events from all channels
    tic3 = tic;
    M = zeros(4,length(timestamps{ch}));
    for i=active_ch_for_tt
        M(i,spikeidx_all_ch{i})=csc{i,1}(spikeidx_all_ch{i});
    end
    M = max(abs(M));
    
    [~,spikeidx] = findpeaks(M,'MinPeakHeight',min(thresh(:)));
    spiketime = timestamps{ch}(spikeidx);
    clear M;
    clear xx;
    clear idx;
    disp(['2) remove redundent spike events on TT',num2str(ii_tetrode),' elepsed:',num2str(toc(tic3)/60),' min']);
    
    %% calculate spikewave
   % spikeidx=[spikeidx_all_ch{:}];
    tic4=tic;
    spikewave = zeros(32,size(spikeidx,2),4);
    mat= repmat([-7:24]',1,length(spikeidx));
    spike_idx = repmat(spikeidx,[32,1]);
    spike_idx = spike_idx+mat;
    
    for ch=active_ch_for_tt
        spikewave(:,:,ch) = csc{ch}(spike_idx);
    end
    
    disp(['3) calculating waveforms TT',num2str(ii_tetrode),' (',num2str(size(spiketime,2)),' spikes)',' elepsed:',num2str(toc(tic4)/60),' min']);
    %clear csc;
    
    if create_movie
        video_file_name=[base_output_dir,'/threshold_mov_TT',num2str(ii_tetrode)];
        create_thresh_movie(video_file_name,timestamps,spiketime,spikewave,thresh)
    end
    
    %% compare to acceptable spike shapes
    if compare_to_library==1
       tic5 = tic;

        [~,iimax ] = max(squeeze(abs(spikewave(8,:,:))),[],2);
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
                
    end
   disp(['4) Compare to library TT',num2str(ii_tetrode),' (',num2str(size(waves,2)),' spikes)',' elepsed:',num2str(toc(tic5)/60),' min']);

   %% Remove spikes within a window (avoid duplications)
   tic6=tic;
    M = zeros(4,length(timestamps{ch}));
    for i=active_ch_for_tt
        M(i,spikeidx(acceptedidx))=csc{i,1}(spikeidx(acceptedidx));
    end
    M = max(abs(M));
    
    [~,spikeidx] = findpeaks(M,'MinPeakDistance',min_sep_events,'MinPeakHeight',min(thresh(:)));
    spiketime = timestamps{ch}(spikeidx);
    clear M;
    clear xx;
    clear idx;
    %% calculate spikewave
    spikewave = zeros(32,size(spiketime,2),4);
    mat= repmat([-7:24]',1,length(spikeidx));
    spike_idx = repmat(spikeidx,[32,1]);
    spike_idx = spike_idx+mat;
    
    for ch=active_ch_for_tt
        spikewave(:,:,ch) = csc{ch}(spike_idx);
    end
    
    
    Timestamps_accepted_spikes_TT{ii_tetrode}=spiketime(:)';
    spikes_TT{ii_tetrode}=spikewave(:,:,:);
      
    
   disp(['5) Remove spikes within a time window (avoid duplications) TT',num2str(ii_tetrode),' (',num2str(length(spike_idx)),' spikes)',' elepsed:',num2str(toc(tic6)/60),' min']);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Coincidence-Detection across Tetrodes to eliminate artefacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic7 = tic;
if length(use_for_sorting)>1
    idx_coincidence_vec = cell(length(use_for_sorting),1);
    
    for ii_spikes = 1:length(Timestamps_accepted_spikes_TT{1}) % Loop over the spikes of the FIRST tetrode
        temp_stack = cell(length(use_for_sorting),1);
        t_spike = Timestamps_accepted_spikes_TT{1}(ii_spikes);
        temp_stack{1,1} = ii_spikes;
        for TT=2:length(use_for_sorting)
            temp_stack{TT,1} = find( abs( Timestamps_accepted_spikes_TT{TT} - t_spike ) <= coincidence_window ); % THE COINCIDENCE DETECTION
        end
        ii = ~cellfun(@isempty,temp_stack);
        if sum(ii)>=length(use_for_sorting)
            for jj=1:length(use_for_sorting)
                idx_coincidence_vec{jj,1} = cat(2,idx_coincidence_vec{jj,1}, temp_stack{jj,1});
            end
        end
    end
end

% %
for i=1:length(use_for_sorting)
    idx_coincidence_vec{i} = unique(idx_coincidence_vec{i});
    if ~isempty(idx_coincidence_vec{i})
        Timestamps_accepted_spikes_TT{i}(idx_coincidence_vec{i})=[];
        spikes_TT{i}(:,idx_coincidence_vec{i},:)=[];
    end
end
disp(['6) cross TT coincident spike removal...',' elepsed:',num2str(toc(tic7)/60),' min']);


%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the data in Neuralynx NTT files for three cases:
%%%%%%%%%%%%%%%%%%%%%%%%%

tic8 = tic;
for ii_Tetrode=1: length(use_for_sorting)
    
    %base_output_dir						= [param.path.analysis, 'preprocessing_spikes4',filesep];
    if ~exist(base_output_dir,'dir'), mkdir(base_output_dir); end
    filename_out = [base_output_dir, 'spikes__TT', num2str(use_for_sorting(ii_Tetrode))];
    
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
disp(['7) Saving spikes to disk...',' elepsed:',num2str(toc(tic8)/60),' min']);
disp(['8) Total detection on all TTs elepsed:',num2str(toc(tic1)/60),' min']);





