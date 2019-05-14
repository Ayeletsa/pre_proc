function PRE_create_cells_struct (p_in,rows)

C = PRE_read_excel_sheet(p_in.excel_sheet ,'Cells',rows,p_in.numeric_fields,[]);
P = NLG_PRE_read_and_fill_excel_sheet (p_in,p_in.excel_sheet,'Experiments',rows);

for ii_cell = 1:length(C)
   day_structs_folder=fullfile(P(1).path_dataout,P(1).year_bat_path,'day_structs');

    day_struct_name = sprintf ('bat_%d_day_%d',C(ii_cell).bat,C(ii_cell).day);
    file_name = fullfile(day_structs_folder,day_struct_name);
    p = load (file_name);
    p = p.p;
    
    data_out_dir = p.path_dataout;
    
    %% 1. load spike file:
    filename_spike = fullfile(data_out_dir,P(1).year_bat_path,...
        num2str(C(ii_cell).day),'detected_spikes',['spikes__TT',num2str(C(ii_cell).TT),'_SS_0',num2str(C(ii_cell).cell_id),'.ntt']);
    
    % ==== Read the Ntt files, and extract the data that occureed within the session: ====
    FieldSelection = [1 1 1 1 1] ; % Will read ALL the parameters, including the Samples (waveforms)
    ExtractionMode = 1 ; % Mode 1 = "Extract All"
    ExtractionModeArray = [] ; % Will read all the data
    
    [Timestamps_usec, ChanNum, CellNumbersSpikeSorting, NumValidSamples, Samples, NlxHeader] = ...
        Nlx2MatSpike( filename_spike, FieldSelection, 1, ExtractionMode, ExtractionModeArray ) ;
    
    %% 2. fit timestamps to bsp
    
    %a. for spikes
    Timestamps_msec = Timestamps_usec / 1e3;
    Timestamps_fitted_to_bsp_msec = polyval (p.time_conv_p_msec_nlx2bsp,Timestamps_msec);
   
    %b. for nlg events
    start_time_msec = [p.S.start_time] / 1e3;
    start_time_fitted_to_bsp_msec =num2cell( polyval (p.time_conv_p_msec_nlx2bsp,start_time_msec)');
    [p.S(:).start_time_fitted_to_bsp_msec]= deal(start_time_fitted_to_bsp_msec{:});
    
    end_time_msec =[ p.S.end_time ]/ 1e3;
    end_time_fitted_to_bsp_msec = num2cell(polyval (p.time_conv_p_msec_nlx2bsp,end_time_msec));
    [p.S(:).end_time_fitted_to_bsp_msec]= deal(end_time_fitted_to_bsp_msec{:});
    
    
    %% 3. load TT info and compute cells separation
    
    L_Ratio = 0; Isolation_dis = 0;
    NTT_file_all =   fullfile(data_out_dir,P(1).year_bat_path,...
        num2str(C(ii_cell).day),'detected_spikes',['spikes__TT',num2str(C(ii_cell).TT),'.ntt']);
    
    [all_spike_ts,all_spike_samples] = Nlx2MatSpike(NTT_file_all, [1 0 0 0 1], 0, 1, []);
    
    
    % Maya's way to compute cluster quality (for test)
    % %         ind = [];
    % %         for n=1:length(spikes_ts)
    % %             ind = [ind find(all_spike_ts == spikes_ts(n))];
    % %         end
    % %         Fet = Create_FeatureSpace(permute(all_spike_samples,[3,2,1]));
    % %         [CluSep, m] = Cluster_Quality_maya(Fet, ind);
    
    [L_Ratio,Isolation_dis] = cluster_quality(all_spike_samples,all_spike_ts,Timestamps_usec);
    
    %% 4.add all important data to cell struct
    
    %a. cell information:
    %--------------------------------------------------------------
    % TODO: change names of params if changed in P
    cell_struct.cell_info=C(ii_cell);
    cell_struct.cell_info.TT_depth=p.depth(C(ii_cell).TT);
    %TO DO - take the threshold from spike detection!
    % cell_struct.cell_info.Spike_threshold_uV_units=p.Spike_threshold_uV_units(:,C(ii_cell).TT);
    cell_struct.cell_info.include_negative_threshold=p.include_negative_threshold;
    cell_struct.cell_info.compare_to_library=p.compare_to_library  ;
    cell_struct.cell_info.do_coincidence_detection=p.do_coincidence_detection;
    cell_struct.cell_info.r_threshold=p.r_threshold  ;
    cell_struct.cell_info.threhold_factor=p.threhold_factor; %for voltage threshold (threhold=median*factor)
    cell_struct.cell_info.min_sep_events=p.min_sep_events; %min separation between events - ~=500us
    cell_struct.cell_info.coincidence_window=p.coincidence_window; % Duration of coincidence-detection window (in microsec) ACROSS tetrodes
    cell_struct.cell_info.thresh_for_coinc=p.thresh_for_coinc; %remove if X tetrodes has event in the window
%    if ~isempty(p.voltage_treshold_for_detection)
%     cell_struct.cell_info.voltage_treshold_for_detection=p.voltage_treshold_for_detection(:,C(ii_cell).TT);
%    end
    cell_struct.cell_info.reference_channel=p.reference_channel;
    cell_struct.cell_info.time_conv_p_msec_nlx2bsp=p.time_conv_p_msec_nlx2bsp;
    %b. experiment information:
    %--------------------------------------------------------------
    cell_struct.exp_info.bat=p.bat;
    cell_struct.exp_info.day=p.day;
    cell_struct.exp_info.path_day_dir=p.path_day_dir;
    cell_struct.exp_info.path_datain=p.path_datain;
    cell_struct.exp_info.path_dataout=p.path_dataout;
    cell_struct.exp_info.bsp_tag_self=p.bsp_tag_self;
    cell_struct.exp_info.bsp_tag_other=p.bsp_tag_other;
    %cell_struct.exp_info.landmarks = p.landmarks;
    cell_struct.exp_info.throw_away_times=p.throw_away_times;
    cell_struct.exp_info.nlg_events=p.S;
    %c. bsp data
    %--------------------------------------------------------------
    cell_struct.bsp_data=p.bsp_data;
    
    %d.spikes:
    %--------------------------------------------------------------
    cell_struct.spikes.spikes_ts_usec=Timestamps_usec;
    cell_struct.spikes.spikes_ts_msec=Timestamps_msec;
    cell_struct.spikes.spikes_ts_fitted_to_bsp_msec=Timestamps_fitted_to_bsp_msec;
    cell_struct.spikes.Samples=Samples;
    cell_struct.spikes.L_Ratio=L_Ratio;
    cell_struct.spikes.Isolation_dis=Isolation_dis;
    
       % calculate stabillity
%         min_bin=2;
% 
%     for behav_i=1:length(cell_struct.exp_info.nlg_events)
%         start_ts_msec=cell_struct.exp_info.nlg_events(behav_i).start_time_fitted_to_bsp_msec;
%         end_ts_msec=cell_struct.exp_info.nlg_events(behav_i).end_time_fitted_to_bsp_msec;
%         time_bins_vec{behav_i}=start_ts_msec:min_bin*60*1e3:end_ts_msec;
%         relevant_spikes=Timestamps_fitted_to_bsp_msec(find(Timestamps_fitted_to_bsp_msec> start_ts_msec & Timestamps_fitted_to_bsp_msec<end_ts_msec ));
%         fr_hist{behav_i}=hist(relevant_spikes,time_vec{behav_i})/(min_bin*60); %hz
%     end
%     
%     pre_sleep_fr = mean(fr_hist{1});
%     post_sleep_fr = mean(fr_hist{end});
%     stability_index = max([pre_sleep_fr post_sleep_fr]) / min([pre_sleep_fr post_sleep_fr]);
%     stability_during_sleep = {pre_sleep_fr,post_sleep_fr,stability_index};
%     cell_struct.spikes.stability = {time_bins_vec,fr_hist,stability_during_sleep};
    
    %5. Save cell struct
    %TO DO - change to proc_data dir...
    stuct_file_name=fullfile(data_out_dir,P(1).year_bat_path,'cell_structs',['bat_',num2str(p.bat),'_day_',num2str(p.day),'_cell_',num2str(cell_struct.cell_info.cell_num),'.mat']);
    
    save(stuct_file_name,'cell_struct')
    
    
end


end