function PRE_create_cells_struct (p_in,rows)

C = PRE_read_excel_sheet(excel_sheet,'Cells',rows,p_in.numeric_fields,[]);

for ii_cell = 1:length(C)
    
    date_ind = find ([P.day] == C(ii_cell).day);
    bat_ind = find ([P.bat ]== C(ii_cell).bat);
    exp_ind = intersect(date_ind,bat_ind);
    
    data_out_dir = P(exp_ind).path_dataout;
    
    %% 1. load spike file:
    filename_spike = fullfile(data_out_dir,['BAT',num2str(C(ii_cell).bat)],...
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
    Timestamps_fitted_to_bsp_msec = polyval (P(exp_ind).time_conv_p_msec_nlx2bsp,Timestamps_msec);
    
    %b. for nlg events
    start_time_msec = [P(exp_ind).S.start_time] / 1e3;
    start_time_fitted_to_bsp_msec =num2cell( polyval (P(exp_ind).time_conv_p_msec_nlx2bsp,start_time_msec)');
    [P(exp_ind).S(:).start_time_fitted_to_bsp_msec]= deal(start_time_fitted_to_bsp_msec{:});
    
    end_time_msec =[ P(exp_ind).S.end_time ]/ 1e3;
    end_time_fitted_to_bsp_msec = num2cell(polyval (P(exp_ind).time_conv_p_msec_nlx2bsp,end_time_msec));
    [P(exp_ind).S(:).end_time_fitted_to_bsp_msec]= deal(end_time_fitted_to_bsp_msec{:});
    
    
    %% 3. load TT info and compute cells separation
    
    L_Ratio = 0; Isolation_dis = 0;
    NTT_file_all = fullfile(data_out_dir,['BAT',num2str(C(ii_cell).bat)],...
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
    % cell_struct.cell_info.TT_depth=P(exp_ind).depth(C(ii_cell).TT);
    %TO DO - take the threshold from spike detection!
    % cell_struct.cell_info.Spike_threshold_uV_units=P(exp_ind).Spike_threshold_uV_units(:,C(ii_cell).TT);
    cell_struct.cell_info.include_negative_threshold=P(exp_ind).include_negative_threshold;
    cell_struct.cell_info.do_comparison_acceptable_shapes=P(exp_ind).do_comparison_acceptable_shapes  ;
    cell_struct.cell_info.do_coincidence_detection=P(exp_ind).do_coincidence_detection;
    cell_struct.cell_info.r_threshold=P(exp_ind).r_threshold  ;
    cell_struct.cell_info.reference_channel=P(exp_ind).reference_channel;
    cell_struct.cell_info.time_conv_p_msec_nlx2bsp=P(exp_ind).time_conv_p_msec_nlx2bsp;
    %b. experiment information:
    %--------------------------------------------------------------
    cell_struct.exp_info.bat=P(exp_ind).bat;
    cell_struct.exp_info.day=P(exp_ind).day;
    cell_struct.exp_info.path_day_dir=P(exp_ind).path_day_dir;
    cell_struct.exp_info.path_datain=P(exp_ind).path_datain;
    cell_struct.exp_info.path_dataout=P(exp_ind).path_dataout;
    cell_struct.exp_info.bsp_tag_self=P(exp_ind).bsp_tag_self;
    cell_struct.exp_info.bsp_tag_other=P(exp_ind).bsp_tag_other;
    cell_struct.exp_info.throw_away_times=P(exp_ind).throw_away_times;
    cell_struct.exp_info.nlg_events=P(exp_ind).S;
    %c. bsp data
    %--------------------------------------------------------------
    cell_struct.bsp_data=P(exp_ind).bsp_data;
    
    %d.spikes:
    %--------------------------------------------------------------
    cell_struct.spikes.spikes_ts_usec=Timestamps_usec;
    cell_struct.spikes.spikes_ts_msec=Timestamps_msec;
    cell_struct.spikes.spikes_ts_fitted_to_bsp_msec=Timestamps_fitted_to_bsp_msec;
    cell_struct.spikes.Samples=Samples;
    cell_struct.spikes.L_Ratio=L_Ratio;
    cell_struct.spikes.Isolation_dis=Isolation_dis;
    % calculate stabillity
        min_bin=2;
    bar_colors=[0.5 0.5 0.5; 0 0 1;0 0 1;0.5 0.5 0.5];
    axes('position',[position_x(3) position_y(1) panel_size/2])
    for behav_i=1:length(cell_struct.exp_info.nlg_events)
        start_ts_msec=cell_struct.exp_info.nlg_events(behav_i).start_time_fitted_to_bsp_msec;
        end_ts_msec=cell_struct.exp_info.nlg_events(behav_i).end_time_fitted_to_bsp_msec;
        
        time_vec{behav_i}=start_ts_msec:min_bin*60*1e3:end_ts_msec;
        relevant_spikes=spikes_ts_fitted_to_bsp_msec(find(spikes_ts_fitted_to_bsp_msec> start_ts_msec & spikes_ts_fitted_to_bsp_msec<end_ts_msec ));
        spike_counts{behav_i}=hist(relevant_spikes,time_vec{behav_i})/(min_bin*60); %hz
        
        bar(time_vec{behav_i},h{behav_i},'FaceColor',bar_colors(behav_i,:))
        hold all
        max_fr_behav(behav_i)=max(h{behav_i});
    end
    
    max_fr=max(max_fr_behav);
    plot([obstacle_time obstacle_time],[0 max_fr],'r')
    title('Firing rate over time')
    ylabel('Firing rate (Hz)')
    xlabel('Time (min)')
    set(gca,'Ylim',[0 max_fr*1.1])
    set(gca,'xlim',[min(time_vec{1}) max(time_vec{end})])
    tick=get(gca,'xtick');
    set(gca,'xtick',tick,'xticklabel',round((tick/1e3)/60))

     %e.Audio:
    %--------------------------------------------------------------
   

    %5. Save cell struct
    %TO DO - change to proc_data dir...
    stuct_file_name=fullfile(data_out_dir,['BAT',num2str(C(ii_cell).bat)],'cell_structs',['bat_',num2str(P(exp_ind).bat),'_day_',num2str(P(exp_ind).day),'_cell_',num2str(cell_struct.cell_info.cell_num),'.mat']);

    save(stuct_file_name,'cell_struct')
    
    
end


end