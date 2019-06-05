function PRE_create_cells_struct (p_in)

cell_rows=p_in.cell_rows;

C = PRE_read_excel_sheet(p_in.excel_sheet ,'Cells',cell_rows,p_in.numeric_fields,[]);
P = NLG_PRE_read_and_fill_excel_sheet (p_in,p_in.excel_sheet,'Experiments',cell_rows);

for ii_cell = 1:length(C)
   day_structs_folder=fullfile(P(1).path_dataout,P(1).year_bat_path,'day_structs');

    day_struct_name = sprintf ('bat_%d_day_%d',C(ii_cell).bat,C(ii_cell).day);
    file_name = fullfile(day_structs_folder,day_struct_name);
    p = load (file_name);
    p = p.p;
    
    main_dir = p.path_day_dir;
    date = char (regexp (main_dir,'\d{8}','match'));
    data_out_dir = p.path_dataout;
    main_dir_out = [p.path_dataout,'\',p.year_bat_path,'\',date];
    dir_sorting = fullfile(main_dir_out,'\spike_sorting\');

    %% 1. load spike file:
    %% Get files
    TT=C(ii_cell).TT;
    spike_TT_folder = dir( fullfile(dir_sorting,['*_TT',num2str(TT),'*']) );
    spike_TT_file=fullfile(dir_sorting,spike_TT_folder.name);
    cell_id=C.cell_id;
    
    % ==== Read the Ntt files, and extract the data that occureed within the session: ====
    FieldSelection = [1 1 1 1 1] ; % Will read ALL the parameters, including the Samples (waveforms)
    ExtractionMode = 1 ; % Mode 1 = "Extract All"
    ExtractionModeArray = [] ; % Will read all the data
    
    [Timestamps_usec, ChanNum, CellNumbersSpikeSorting, NumValidSamples, Samples, NlxHeader] = ...
        Nlx2MatSpike( spike_TT_file, FieldSelection, 1, ExtractionMode, ExtractionModeArray ) ;
    %Take only relevant spikes for cell
    Timestamps_usec=Timestamps_usec(CellNumbersSpikeSorting==cell_id);
    Samples=Samples(:,:,CellNumbersSpikeSorting==cell_id);
    
    %% 2. fit timestamps to bsp
    
%     %a. for spikes
%     Timestamps_msec = Timestamps_usec / 1e3;
%     Timestamps_fitted_to_bsp_msec = polyval (p.time_conv_p_msec_nlx2bsp,Timestamps_msec);
%    
%     %b. for nlg events
%     start_time_msec = [p.S.start_time] / 1e3;
%     start_time_fitted_to_bsp_msec =num2cell( polyval (p.time_conv_p_msec_nlx2bsp,start_time_msec)');
%     [p.S(:).start_time_fitted_to_bsp_msec]= deal(start_time_fitted_to_bsp_msec{:});
%     
%     end_time_msec =[ p.S.end_time ]/ 1e3;
%     end_time_fitted_to_bsp_msec = num2cell(polyval (p.time_conv_p_msec_nlx2bsp,end_time_msec));
%     [p.S(:).end_time_fitted_to_bsp_msec]= deal(end_time_fitted_to_bsp_msec{:});
     
    
    %% 3. load TT info and compute cells separation
    
    L_Ratio = 0; Isolation_dis = 0;
  
    
    [all_spike_ts,all_spike_samples] = Nlx2MatSpike(spike_TT_file, [1 0 0 0 1], 0, 1, []);
    
    
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
     
    %b. experiment information:
    %--------------------------------------------------------------
    cell_struct.exp_info.bat=p.bat;
    cell_struct.exp_info.day=p.day;
    cell_struct.exp_info.path_day_dir=p.path_day_dir;
    cell_struct.exp_info.path_datain=p.path_datain;
    cell_struct.exp_info.path_dataout=p.path_dataout;
    cell_struct.exp_info.bsp_tag_self=p.bsp_tag_self;
    cell_struct.exp_info.bsp_tag_other=p.bsp_tag_other;
    cell_struct.exp_info.landmarks = p.landmarks;
    cell_struct.exp_info.throw_away_times=p.throw_away_times;
    cell_struct.exp_info.nlg_events=p.S;
    
    %c. bsp data
    %--------------------------------------------------------------
    cell_struct.bsp_data=p.bsp_data;
    
    %d.spikes:
    %--------------------------------------------------------------
    cell_struct.spikes.spikes_ts_usec=Timestamps_usec;
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
    stuct_file_name=fullfile(P.path_dataout,P.year_bat_path  ,'cell_structs',['bat_',num2str(p.bat),'_day_',num2str(p.day),'_cell_',num2str(cell_struct.cell_info.cell_num),'.mat']);
    
    save(stuct_file_name,'cell_struct')
    
    
end


end