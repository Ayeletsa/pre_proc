function p = BSP_PRE_extract_data (p)

main_dir = p.path_day_dir;
bsp_dir = fullfile(main_dir,'bsp');

% extract bsp data from raw recordings data
bsp_data = bsp_extract_data (bsp_dir);
% use only wanted tags data and delete the rest
self_other_bsp_tags = [p.bsp_tag_self p.bsp_tag_other];
all_recorded_tags = [bsp_data.tag_ID];
bsp_tags_ind = ismember(all_recorded_tags,self_other_bsp_tags);
bsp_data = bsp_data(bsp_tags_ind);
%%%%
% calib=load('calib_tunnel.mat');
% calib=calib.calib_tunnel;
% pos_linearized = POS_calc_linearized([bsp_data(1).pos(:,1:2)],calib);
%%%%%
% normelize the data (bring to origin, rotate to correct axis)
[bsp_data] = normelize_my_data2 (bsp_data);
% clean unwanted samples and smooth the data
bsp_data = clean_and_smooth_bsp_data(bsp_data);
% interpolate locaion to a common time vector
bsp_data = upsample_100hz_modified(bsp_data);

p.bsp_data = bsp_data;
%p.landmarks = landmarks;

end