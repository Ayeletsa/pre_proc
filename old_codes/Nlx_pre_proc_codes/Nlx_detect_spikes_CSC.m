function Nlx_detect_spikes_CSC(expname,forcerecalc)

%% IMPORT
%
% Didi, 
% November 2015 
% Adopted from Arseny and Michael 
% 
% Here we detect the spikes for the single tetrode on each antenna and save them.
% Output is 2 NTT files containing:
% 1. detected spikes
% 2. spikes that were detected but thrown away after library comparison

eval(expname);

folder_name = param.path.spikes_detection;
if exist(folder_name,'dir') && ~forcerecalc
    disp(['Detect-spikes was already done',expname]);
    return; 
end

%% detect spikes
for ii_tetrode=1:length(param.tetrodes.use_tetrodes)
    TT = param.tetrodes.use_tetrodes(ii_tetrode)
    %% load raw data (CSCs)
    csc = {};
    timestamps  = {};
    for i_ch=1:4
        filename =  fullfile(param.path.CSC_spikes,['spikes_b',param.bat,'_d',param.day,'_TT',num2str(TT),'_ch',num2str(i_ch-1),'.ncs']);
        [csc{i_ch}, timestamps{i_ch}] = Nlx_csc_read(filename,param.spikes.neural_data_ts);
    end
    % sanity check - make sure all csc files are in the same length
    if range([cellfun(@length, timestamps)]) == 0
        timestamps = timestamps{1};
    else
        error('CSC file of channels from the same TT are different in length!!!')
    end
    
        
    %% detect threshold crossing 
    disp('------------------------------------')
    disp('Detect spike (thr crossing)         ')
    tic
    thres_cross_vec ={};
    thres_cross_IX  ={};
    SPK_start = {}; 
    SPK_End = {}; 
    
    Last_Spike_IX = 0; % Initialize
    SPK_timestamp = [];% Initialize

    % TODO: run only over valid channels
    for i_ch = 1:4   
        spike_thr_uV_units = param.spikes.spike_threshold_uV_units(TT,i_ch);
        thres_cross_vec{i_ch} = zeros(1,length(csc{i_ch}));
        if param.spikes.include_negative_threshold
            thres_cross_IX{i_ch} = find(csc{i_ch}>spike_thr_uV_units | csc{i_ch}<-spike_thr_uV_units);
        else
            thres_cross_IX{i_ch} = find(csc{i_ch}>spike_thr_uV_units);
        end;
        
        % place '1' at threshold crossing:
        thres_cross_vec{i_ch}(thres_cross_IX{i_ch})=1; % place '1' at threshold crossing

        % Find the start and end segment of each spiking event on each channel of the tetrode:
        diff_thres_cross_vec{i_ch} = diff(thres_cross_vec{i_ch});
        SPK_start{i_ch} = find(diff_thres_cross_vec{i_ch} ==1)+1;
        SPK_End{i_ch} = find(diff_thres_cross_vec{i_ch} ==(-1))+1;
        
        %Check for the case that csc start/end point already crossed the thr
        if SPK_End{i_ch}(1) < SPK_start{i_ch}(1)
            SPK_End{i_ch}(1) = [];
        end
        if SPK_start{i_ch}(end) > SPK_End{i_ch}(end)
            SPK_start{i_ch}(end) = [];
        end

        % correct the spikes index according to maximum
        temp_SPK_events = []; % Initialize
        SPK_max{i_ch} = [];
        SPK_IX{i_ch} = []; % with regard to the maximum value of events that crossed thr
        for jj = 1:length(SPK_End{i_ch})
            temp_IX_vec = SPK_start{i_ch}(jj):SPK_End{i_ch}(jj);
            temp_SPK_events = csc{i_ch}(temp_IX_vec);
            [temp_SPK_max,IX] = max(temp_SPK_events);
            max_IX = temp_IX_vec(IX);
            SPK_IX{i_ch}(1,jj) = max_IX; % Correct for the real timestamp relative to the entire recording session
        end
        
    end
    clear diff_thres_cross_vec
    clear thres_cross_vec
    clear thres_cross_IX
    toc
  
    
    %% Merge spikes from different channels
    % Merge spikes from all channels of the same TT
    % If two spikes are detected on differnt channel but they are close
    % enough in time (<x_sep_spike_thres), we merge them to a single
    % spike
    disp('------------------------------------')
    disp('Merge spikes from different channels')
    tic
    
% % % % % % % % % % % % % % % % % % % % % % % %     
    events_all_ch = zeros(4,length(csc{1}));
    for ii_ch = 1:4
        events_all_ch( ii_ch, SPK_IX{ii_ch} ) = 1;
    end
    events_ch_combined = sum(events_all_ch);
    clear events_all_ch;
    toc
% % % % % % % % % % % % % % % % % % % % % % % %     
    rect_kernel = ones(1,param.spikes.x_sep_spike_thres);
    events_ch_combined_filtered = filtfilt(rect_kernel, 1, events_ch_combined);
    toc
% % % % % % % % % % % % % % % % % % % % % % % %     
    [~,SPK_IXs_All_Ch_combined_sorted] = findpeaks(events_ch_combined_filtered);
%     [~,SPK_IXs_All_Ch_combined_sorted] = findpeaks(events_ch_combined, 'MinPeakDistance', param.spikes.x_sep_spike_thres);
    clear events_ch_combined
    toc
% % % % % % % % % % % % % % % % % % % % % % % %     


% % %     % take as a start point all the spikes from ch1
% % %     SPK_All_Ch_combined = SPK_IX{1};
% % % 
% % %     % merge into this list spikes from all other channels
% % %     for curr_chan_idx=2:4
% % %         % First loop over the spikes detected on the first channel and compare their seperation from those detected on the second channel
% % %         % Save only those which are seperated by a minimal number of bins to avoid counting the same spike twice
% % %         tic
% % %         for ii_spike = 1:length(SPK_All_Ch_combined)
% % %             current_spike_IX = SPK_All_Ch_combined(ii_spike);
% % %             shared_IXs = find(abs(SPK_IX{curr_chan_idx} - current_spike_IX)<= param.spikes.x_sep_spike_thres);
% % %             if ~isempty(shared_IXs) % i.e., the same spike is detected twice
% % %                 SPK_IX{curr_chan_idx}(shared_IXs) = [];
% % %             end
% % %         end
% % %         toc
% % %         % Now that we do not have the same spikes on two channels we can merge the spike IXs of both channels,
% % %         % as those represent different spikes
% % %         SPK_All_Ch_combined = [SPK_All_Ch_combined,SPK_IX{curr_chan_idx}];
% % %     end
% % %     
% % %     % sort them to preserve their temporal order:
% % %     SPK_IXs_All_Ch_combined_sorted = sort(SPK_All_Ch_combined);
    
    toc
    
    %% Extract spike waveforms
    % Extract the waveform and timestamps of each spike (as in Neuralynx, the peak
    % will be the 8th sample out of over all 32 samples of
    %the spike shape vec.
    
    disp('------------------------------------')
    disp('Extract waveforms')
    tic
    
    SPK_waveforms = zeros(4,length(SPK_IXs_All_Ch_combined_sorted),32);% Initialize
    current_file_spike_counter = 0;
    for ii_spike = 1:length(SPK_IXs_All_Ch_combined_sorted)
        current_spike_max_IX = SPK_IXs_All_Ch_combined_sorted(ii_spike);
        if ((current_spike_max_IX+24<=length(timestamps))&& (current_spike_max_IX-7>0)) % i.e., we are NOT cutting the spike in the middle, TODO: no need to check this every loop, we can check only for the first and last spikes....
            current_file_spike_counter = current_file_spike_counter + 1;
            SPK_timestamp(1,Last_Spike_IX+current_file_spike_counter) = timestamps(current_spike_max_IX);
            for curr_chan_idx=1:4
                SPK_waveforms(curr_chan_idx,Last_Spike_IX+current_file_spike_counter,:) = csc{curr_chan_idx}(current_spike_max_IX-7:1:current_spike_max_IX+24);
            end
        end
    end
    Last_Spike_IX = Last_Spike_IX + current_file_spike_counter;
    
    toc
    
    %% library of acceptable spike shapes
    % Clean artifacts using the library of acceptable spike shapes.
    % Now we will throw away noisy theshold crossing using the library of
    % acceptable spike shapes as reference

    disp('------------------------------------')
    disp('library of acceptable spike shapes')
    tic
    
    % First we will need to normalize peak amplitude of each spike to a value
    % of '1' such that we can comapre it to the library of acceptable spike
    % shapes:
    
    load(param.path.spike_shapes_lib);
    vector_of_accepted_spikes = zeros( 1, length(SPK_waveforms) ) + NaN ; % Initialize
    vector_of_max_r_values = zeros( 1, 1 ) + NaN ;
    r_threshold = param.spikes.r_threshold;
    
    % For each event take the channel with the largest peak (8th point)
    spikes_waveforms = zeros(size(SPK_waveforms,2),32);
    for ii_spike = 1:size(SPK_waveforms,2)
        [~,ii_ch_max] = max(SPK_waveforms(:,ii_spike,8));
        spikes_waveforms(ii_spike,:) = SPK_waveforms(ii_ch_max,ii_spike,:);
    end
    
    % calc corr
    xxx_lags_shifts = [1:30; 2:31; 3:32];
    ccc = [];
    rrr = [];
    for ii_shift = 1:size(xxx_lags_shifts,1)
        xxx_lags = xxx_lags_shifts(ii_shift, :);
        ccc = corr(spikes_waveforms(:,xxx_lags)', library_of_acceptable_spike_shapes(:,xxx_lags)');
        rrr(ii_shift,:) = max(ccc,[],2);
    end
    rrr = max(rrr,[],1);
    
    vector_of_accepted_spikes = ( rrr >=  r_threshold ); % TODO: plot hist of 'r' values + thr
    
% % % %     for ii_spike = 1:size(SPK_waveforms,2)
% % % %         %(ii_spike/size(SPK_waveforms,2))*100
% % % %         spike_shape_4channels = zeros(32,4);
% % % %         for jj_channel = 1:size(SPK_waveforms,1)
% % % %             spike_shape_4channels(:,jj_channel) = SPK_waveforms{jj_channel,ii_spike}';
% % % %         end
% % % %         % Choose the channel # for which the spike has the largest height:
% % % %         [ stam  idx_channel_max_height ] = max( max( spike_shape_4channels ) );
% % % %         spike_shape = spike_shape_4channels( :, idx_channel_max_height )' ;
% % % %         
% % % %         if ( std( spike_shape(2:end-1) ) == 0 ), % If this is a completely FLAT "spike", I cannot compute CORRCOEF, so I will set r = 0
% % % %             vector_of_max_r_values( ii_spike ) = 0 ;  % Set r = 0 in this case
% % % %         else % If this spike DOES have some shape (this is the case basically for ALL the recorded waveforms)
% % % %             % Compute the correlation coefficients with all the acceptable spike shapes -- lag 0 :
% % % %             xxx_lags = 2 : 31 ; % Use CENTRAL 30 points
% % % %             ccc = corrcoef([ spike_shape(2:end-1) ; library_of_acceptable_spike_shapes(:,xxx_lags)]' ); % Correlation coefficients
% % % %             rrr_vec_lag_0 = ccc( 1, 2:end ); % All the correlation coefficients with the acceptable spike shapes
% % % %             % Compute the correlation coefficients with all the acceptable spike shapes -- lag (+1) :
% % % %             xxx_lags = 1 : 30 ; % Use FIRST 30 points (RIGHT shift of the "Acceptable Spike Shapes matrix")
% % % %             ccc = corrcoef([ spike_shape(2:end-1) ; library_of_acceptable_spike_shapes(:,xxx_lags)]' ); % Correlation coefficients
% % % %             rrr_vec_lag_plus1 = ccc( 1, 2:end ); % All the correlation coefficients with the acceptable spike shapes
% % % %             % Compute the correlation coefficients with all the acceptable spike shapes -- lag (-1) :
% % % %             xxx_lags = 3 : 32 ; % Use LAST 30 points (LEFT shift of the "Acceptable Spike Shapes matrix")
% % % %             ccc = corrcoef([ spike_shape(2:end-1) ; library_of_acceptable_spike_shapes(:,xxx_lags)]' ); % Correlation coefficients
% % % %             rrr_vec_lag_minus1 = ccc( 1, 2:end ); % All the correlation coefficients with the acceptable spike shapes
% % % %             % Save the MAXIMAL r value -- use the maximal correlation among the 3 lags (-1,0,+1):
% % % %             vector_of_max_r_values( ii_spike ) = max( [ rrr_vec_lag_0  rrr_vec_lag_plus1  rrr_vec_lag_minus1 ] );
% % % %         end
% % % %         
% % % %         % Determine if this spike should be Accepted (set value to '1') or Rejected (set value to '0'):
% % % %         vector_of_accepted_spikes( ii_spike ) = ...
% % % %             vector_of_max_r_values( ii_spike )  >=  r_threshold ;
% % % %         % Accept the spike shape ('1') if its correlation with ANY of the acceptable shapes
% % % %         % is >= r_threshold ; else, reject the spike ('0').
% % % %         % FOR DBG: plot waveform in 'r' if rejected and 'b' if accepted
% % % %         %     if (vector_of_sccepted_spikes( ii_spike )== 1)
% % % %         %     plot(spike_shape)
% % % %         %     else
% % % %         %     plot(spike_shape,'r')
% % % %         %     end
% % % %         %     title(['spike num = ',num2str(ii_spike), ' ,r val = ' num2str(vector_of_max_r_values( ii_spike ))])
% % % %         %     waitforbuttonpress
% % % %     end % End "Loop over spikes extracted from the Ntt file"
    

    % Find the IXs of the accepted spike waveforms and store the accepted and
    % not-accepted waveforms speratly.
    IX_accepted = find(vector_of_accepted_spikes ==1);
    IX_NO_accepted = find(vector_of_accepted_spikes ==0);
    
    
    toc
    
    
    
    %%
    % Extract the accepted (and not accepted) spikes and define new varialbes:
    Spike_waveforms_accepted = zeros(4,length(IX_accepted),32);
    Spike_waveforms_NO_accepted = zeros(4,length(IX_NO_accepted),32);
    
    counter_accepted = 0;
    counter_NO_accepted = 0;
    for ii = 1:length(vector_of_accepted_spikes)
        if vector_of_accepted_spikes(ii) == 1;
            counter_accepted = counter_accepted + 1;
            
            for curr_chan_idx=1:4
                Spike_waveforms_accepted(curr_chan_idx,counter_accepted,:) = SPK_waveforms(curr_chan_idx,ii,:);
            end
        else
            counter_NO_accepted = counter_NO_accepted + 1;
            for curr_chan_idx=1:4
                Spike_waveforms_NO_accepted(curr_chan_idx,counter_NO_accepted,:) = SPK_waveforms(curr_chan_idx,ii,:);
            end
        end
    end
    
    
    %%
    Timestamps_accepted_spikes = SPK_timestamp(IX_accepted);
    %rotate the timestamps to be 1 x num_records
    %timestamps = rot90(Timestamps_accepted_spikes);
    
    % Clear un-needed variables
    clear IX_NO_accepted IX_NO_accepted_sleep IX_accepted_sleep
    clear SPK_All_Ch_combined
    clear Spike_waveforms_NO_accepted
    clear Spike_waveforms_NO_accepted_sleep Spike_waveforms_accepted_sleep Timestamps_accepted_spikes_sleep
    clear VT VT_Parameters ccc spikes_sleep
    clear vector_of_accepted_spikes vector_of_accepted_spikes_sleep
    clear csc
    clear timestamps 

    [numCh, numRec ~] = size(Spike_waveforms_accepted);
    spikes = zeros(32,4,numRec);
    for rec=1:numRec
        %rec/numRec
        for channel = 1:numCh
            current_channel_waveform = squeeze(Spike_waveforms_accepted(channel,rec,:));
            %for point=1:32
            %             spikes(point, channel, rec) = current_channel_waveform(point)*amplitude_factor;
            %end
            spikes(:, channel, rec) = current_channel_waveform*param.spikes.amplitude_factor;
        end
    end
    
    Timestamps_accepted_spikes_TT{ii_tetrode}=Timestamps_accepted_spikes;
    spikes_TT{ii_tetrode}=spikes;
end %end looping over Tetrodes




%%  Read Timestamp data and use Coincidence-Detection across Tetrodes to eliminate artifacts (just as we do in the wired case): --------
if param.spikes.do_coincidence_detection
    
    disp('-------------------------------------------')
    disp('Coincidence-Detection (eliminate artifacts)')
    tic
    
    idx_coincidence_vec{length(param.tetrodes.use_tetrodes)} = [];  % Initialize this variable (for later)
    
    % Find coincidence-detection events = coincidence-detection on a millisecond-scale:
    for ii_spikes = 1:length(Timestamps_accepted_spikes_TT{1}), % Loop over the spikes of the FIRST tetrode
        t_spike = Timestamps_accepted_spikes_TT{1}(ii_spikes);
        idx{1} = ii_spikes; % Index of this spike
        
        for ii_file = 2:length(param.tetrodes.use_tetrodes) % Loop over the other tetrodes, finding the coincidnce-detection
            idx{ii_file} = find( abs( Timestamps_accepted_spikes_TT{ii_file} - t_spike ) <= param.spikes.coincidence_window ); % THE COINCIDENCE DETECTION
        end
        
        % Check that the coincidence-detection occurred on ALL tetrodes:
        test_variable = 1;
        for ii_file =1: length(param.tetrodes.use_tetrodes) % Loop over all tetrode files
            if ( isempty( idx{ii_file} ) ), % If there is NO coincidence-detection
                test_variable = 0 ;
            end
        end
        
        if length(param.tetrodes.use_tetrodes)>=3 %If there are at least 3 tetrodes
            % Check that the temporal separation between the OTHER two tetrodes also meets the coincidence-detection criterion
            % (this is needed since the spikes on the OTHER two tetrodes may be up to twice-the-time apart, in principle!!! ):
            if ( test_variable == 1 ), % If we passed the first test
                if ( abs( Timestamps_accepted_spikes_TT{2}(idx{2}(1)) - Timestamps_accepted_spikes_TT{3}(idx{3}(1)) ) > param.spikes.coincidence_window  ),
                    test_variable = 0 ; % Reset test_variable if the time between spikes on the OTHER two tetrodes is too large
                end
            end
        end;
        
        % Save the indexes of the coincidence-detection "spikes" to be REMOVED:
        if ( test_variable == 1 ), % If all indexes exist = there IS a coincidence detection on ALL tetrodes
            for ii_file =1:length( param.tetrodes.use_tetrodes) % Loop over all tetrode files
                idx_coincidence_vec{ii_file} = [ idx_coincidence_vec{ii_file}  idx{ii_file} ];
            end
        end
%         if mod(ii_spikes,1000) == 0
%             disp(['Processing coindicedence detection across TT- ' num2str((ii_spikes/length(Timestamps_accepted_spikes_TT{tetrodes(1)}))*100)])
%         end
    end %end looping over spikes for coincidence detection
    
    toc
    
else
    idx_coincidence_vec{1}=[];idx_coincidence_vec{2}=[];idx_coincidence_vec{3}=[];idx_coincidence_vec{4}=[];
end




%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the data in Neuralynx NTT files for three cases:
% (1) All the data.


for ii_Tetrode=1: length(param.tetrodes.use_tetrodes)

    base_output_dir = param.path.spikes_NTT;
    if ~exist(base_output_dir,'dir'), mkdir(base_output_dir); end 
    filename_out = fullfile(base_output_dir, ['spikes_',param.expname,'_TT', num2str(ii_Tetrode)]);
   
    % load generic header to add to NTT files
    load(param.path.NTT_header);
    
    Timestamps_accepted_spikes=Timestamps_accepted_spikes_TT{ii_Tetrode};
%     Timestamps_accepted_spikes(idx_coincidence_vec{ii_Tetrode})=[];
    spikes=spikes_TT{ii_Tetrode};
%     spikes(:,:,idx_coincidence_vec{ii_Tetrode})=[];

    %only export timestamps and data points - Full session
    FieldSelection = [1 0 0 0 1 1];        
    
    Mat2NlxSpike([filename_out,'.NTT'],0,1,[],FieldSelection,Timestamps_accepted_spikes, spikes,Header);  
    
    
end % end looping over tetrodes

