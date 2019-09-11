function p = BSP_PRE_extract_data (p)

%%% parameters - should be outside!!
%%%%% TODO%%%%%%%%%%%
outliers_speed=20; %m/s
outliers_pos_dist_from_midline=2; %m
resample_fs = 100;
clib_file_name='calib_tunnel.mat'; %write the full path!\

ns_factor=1e9;
%%%
main_dir = p.path_day_dir;
bsp_dir = fullfile(main_dir,'bsp');

% extract bsp data from raw recordings data
bsp_data = bsp_extract_data (bsp_dir);
% use only wanted tags data and delete the rest because sometimes due to
% bug there are non relevant tags
self_other_bsp_tags = [p.bsp_tag_self p.bsp_tag_other];
all_recorded_tags = [bsp_data.tag_ID];
bsp_tags_ind = ismember(all_recorded_tags,self_other_bsp_tags);
bsp_data = bsp_data(bsp_tags_ind);
%%%%
for tag_i=1:length(bsp_data)
%% Linearize (project to the tunnel midline curve)
calib=load(clib_file_name);  
calib=calib.calib_tunnel;
bsp_data(tag_i).pos_linearized = POS_calc_linearized([bsp_data(tag_i).pos(:,1:2)],calib);
bsp_data(tag_i).speed = ns_factor.*[0 ;sqrt(sum(diff(bsp_data(tag_i).pos,1,1).^2,2)) ./ diff(bsp_data(tag_i).ts_ns)];

%% identify and remove outliers
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
bsp_data(tag_i).ts_ns_no_outliers=bsp_data(tag_i).ts_ns;
bsp_data(tag_i).ts_ns_no_outliers(bsp_data(tag_i).outliers_IX) = [];

%convert to us:
bsp_data(tag_i).ts_us_no_outliers=bsp_data(tag_i).ts_ns_no_outliers/1e3;
%% Fill holes (interp/exterp)
[bsp_data(tag_i).pos_fill_holes, bsp_data(tag_i).ts_us_fill_holes] = POS_fill_holes_also_in_y(bsp_data(tag_i).pos_linearized_no_outliers', bsp_data(tag_i).ts_us_no_outliers');

%% UP-sample
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


p.bsp_data = bsp_data;
%p.landmarks = landmarks;

end
end