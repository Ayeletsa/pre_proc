function p = BSP_PRE_PROC_data (p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:
% This code analyze the bsp raw data and process it: remove outliers, fill
% holes and upsample the data.

% This code is very similar to Tamir and Shir's code with 2 matin changes:
%1. works on 2 tags (or more)
%2. analyze the data in y as well
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load params file
load(p.bsp_params_file_name);


us_factor=1e6;
%%%
main_dir = p.path_day_dir;
bsp_dir = fullfile(main_dir,'bsp');

%% 1. extract bsp data from raw recordings data
bsp_data = bsp_extract_data (bsp_dir);
% use only wanted tags data and delete the rest because sometimes due to
% bug there are non relevant tags
self_other_bsp_tags = [p.bsp_tag_self p.bsp_tag_other];
all_recorded_tags = [bsp_data.tag_ID];
bsp_tags_ind = ismember(all_recorded_tags,self_other_bsp_tags);
bsp_data = bsp_data(bsp_tags_ind);
%%%%

for tag_i=1:length(bsp_data)
    % sync ts if needed:
    if strcmp(p.sync_to,'bsp')
        bsp_data(tag_i).ts_nlg_usec=bsp_data(tag_i).ts_ns*10^-3; 
    else
        %DEBUG CHECK TS units!
        bsp_data(tag_i).ts_nlg_usec=interp1(p.sync.bsp_ts_for_sync_with_nlg, p.sync.nlg_ts_for_sync_with_bsp, bsp_data(tag_i).ts_ns, 'linear','extrap');
       
    end
    %remove dots that are on the same x value
    ind_rep=find(diff(bsp_data(tag_i).pos(:,1))==0)+1;
    bsp_data(tag_i).pos(ind_rep,:)=[];
    bsp_data(tag_i).ts_ns(ind_rep)=[];
    bsp_data(tag_i).ts_nlg_usec(ind_rep)=[];
    %% 2. Linearize (project to the tunnel midline curve)
    calib=load(clib_file_name);
    calib=calib.calib_tunnel;
    bsp_data(tag_i).pos_linearized = POS_calc_linearized([bsp_data(tag_i).pos(:,1:2)],calib);
    bsp_data(tag_i).speed = us_factor.*[0 ;sqrt(sum(diff(bsp_data(tag_i).pos,1,1).^2,2)) ./ diff(bsp_data(tag_i).ts_nlg_usec)];
    
    %% 3. identify and remove outliers
    % by distance from the tunnel
    bsp_data(tag_i).outliers_pos_IX = find( abs(bsp_data(tag_i).pos_linearized(:,2)) > outliers_pos_dist_from_midline );
    % by high momentray velocity
    bsp_data(tag_i).outliers_speed_IX = find( bsp_data(tag_i).speed > outliers_speed);
    bsp_data(tag_i).outliers_IX = union(bsp_data(tag_i).outliers_pos_IX, bsp_data(tag_i).outliers_speed_IX);
    bsp_data(tag_i).outliers_pos_dist_from_midline = outliers_pos_dist_from_midline;
    bsp_data(tag_i).outliers_speed = outliers_speed;
    
    %remove outliers:
    bsp_data(tag_i).pos_linearized_no_outliers = bsp_data(tag_i).pos_linearized;
    bsp_data(tag_i).pos_linearized_no_outliers (bsp_data(tag_i).outliers_IX,:) = [];
    bsp_data(tag_i).ts_us_no_outliers=bsp_data(tag_i).ts_nlg_usec;
    bsp_data(tag_i).ts_us_no_outliers(bsp_data(tag_i).outliers_IX) = [];
    
    
    %% 4. Fill holes (interp/exterp)
    [bsp_data(tag_i).pos_fill_holes, bsp_data(tag_i).ts_us_fill_holes] = POS_fill_holes_also_in_y(bsp_data(tag_i).pos_linearized_no_outliers', bsp_data(tag_i).ts_us_no_outliers',bsp_dir);
    
    %% 5. UP-sample
    Ts = 1/resample_fs;
    ts_us_upsampled = bsp_data(tag_i).ts_us_fill_holes(1) : Ts*1e6 : bsp_data(tag_i).ts_us_fill_holes(end);
    for dim_i=1:2 %also y
        pos_upsampled{dim_i} = interp1(bsp_data(tag_i).ts_us_fill_holes, bsp_data(tag_i).pos_fill_holes(dim_i,:), ts_us_upsampled);
    end
    pos_upsampled_mat=cell2mat(pos_upsampled);
    pos_upsampled_mat=reshape(pos_upsampled_mat,[],2)';
    
    bsp_data(tag_i).pos_upsampled=pos_upsampled_mat;
    bsp_data(tag_i).ts_us_upsampled = ts_us_upsampled;
    bsp_data(tag_i).fs_upsampled = resample_fs;
    
    %% 6. load landmakrs
    
    [num,txt,raw] =xlsread(landmark_file);
    [xy,distance,t] = distance2curve(...
        calib.curvexy(1:50:end,:),...        % todo: move the sabsampling to the creation of the calibration...
        [num(:,1) num(:,2)],...
        'linear');
    
    landmarks=t.*calib.tunnel_length;
    
    
end
%% 7. save
p.bsp_data = bsp_data;
p.landmarks = landmarks;
end