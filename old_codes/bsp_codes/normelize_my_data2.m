function [bsp_data] = normelize_my_data2 (bsp_data)

%% create params
% [num,txt,raw] = xlsread('E:\Dan\couples_recording_data\landmarks_locations\tag_1930_landmarks.xlsx');
% landmark_names = txt(2:end,1);

balls_pos=[1372 1493]; % measured balls x pos
median_ball1_pos=zeros(2,2);
median_ball2_pos=zeros(2,2);
n_tags = length(bsp_data);
on_the_balls_ind=cell(1,n_tags);


%% find balls positions by both tags
%TODO: record the balls pos one time for all tags and all experiments.

for ii_tag=1:n_tags
    
    x_pos_diff=[0;diff(bsp_data(ii_tag).pos(:,1))];
    time_diff_sec = [0;diff(bsp_data(ii_tag).ts_ns).*1e-9];
    vel_meter_sec = x_pos_diff./time_diff_sec;
    ball1_ind=(bsp_data(ii_tag).pos(:,1)<(balls_pos(1)+3) & vel_meter_sec<2); % index of measurements on the far ball
    ball2_ind=(bsp_data(ii_tag).pos(:,1)>(balls_pos(2)-3) & vel_meter_sec<2); % index of measurements on the close ball
    median_ball1_pos(ii_tag,:)=[median(bsp_data(ii_tag).pos(ball1_ind,1)),median(bsp_data(ii_tag).pos(ball1_ind,2))]; % meadian x,y of the far ball
    median_ball2_pos(ii_tag,:)=[median(bsp_data(ii_tag).pos(ball2_ind,1)),median(bsp_data(ii_tag).pos(ball2_ind,2))]; % median x,y of the close ball
    on_the_balls_ind{ii_tag}=ball1_ind | ball2_ind;
    
end

balls_pos=zeros(2,2);
balls_pos(1,:)=[mean(median_ball1_pos(:,1)),mean(median_ball1_pos(:,2))];
balls_pos(2,:)=[mean(median_ball2_pos(:,1)),mean(median_ball2_pos(:,2))];

% ball1_ind = find(contains(landmark_names(:,1),'ball1'));
% ball1_pos = num(ball1_ind,1:2);
% ball2_ind = find(contains(landmark_names(:,1),'ball2'));
% ball2_pos = num(ball2_ind,1:2);
% balls_pos = [ball1_pos;ball2_pos]


% subtract the coordinates of the first ball from all data, so the first
% ball will be located at the origin

bsp_data(1).pos(:,1)=bsp_data(1).pos(:,1)-balls_pos(1,1); 
bsp_data(2).pos(:,1)=bsp_data(2).pos(:,1)-balls_pos(1,1);
bsp_data(1).pos(:,2)=bsp_data(1).pos(:,2)-balls_pos(1,2); 
bsp_data(2).pos(:,2)=bsp_data(2).pos(:,2)-balls_pos(1,2);

new_balls_pos=balls_pos-[balls_pos(1,:);balls_pos(1,:)];
ball1_new_pos = new_balls_pos(1,:);
ball2_new_pos = new_balls_pos(2,:);


%% Rotate the hamama so it will lie on the x-axis

% a. calculate hamama slope
slope = (ball2_new_pos(2) - ball1_new_pos(2))/(ball2_new_pos(1) - ball1_new_pos(1));

% b. calculate the slope angle
angle_of_slope_rad = atan(slope);

% c. create a clockwise rotation matrix
a = angle_of_slope_rad;
rotation_matrix_clockwise = [cos(a) sin(a);-sin(a) cos(a)];

% d. rotate bsp coordinates and save tham into the bsp_data struct 
% (z pos untouched)
for ii_tag=1:n_tags
    
    xy_pos = bsp_data(ii_tag).pos(:,1:2)';
    rotated_xy_pos = rotation_matrix_clockwise * xy_pos;
    bsp_data(ii_tag).normalized_pos = [rotated_xy_pos' , bsp_data(ii_tag).pos(:,3)];
    
end

% landmarks_locations(:,1) = num(:,1) - balls_pos(1,1);
% landmarks_locations(:,2) = num(:,2) - balls_pos(1,2);
% new_landmarks_locations = [rotation_matrix_clockwise * landmarks_locations']';
% 
% obs_name = ['obs pos' num2str(obs_num)];
% obs_ind = find(contains(landmark_names(:,1),obs_name));
% obs_location = new_landmarks_locations(obs_ind,1:2);
% 
% landmarks.names = landmark_names;
% landmarks.pos = new_landmarks_locations;
% landmarks.obs = obs_location;

end