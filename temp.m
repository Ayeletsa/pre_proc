bat=2336;
day=20190916;
p.bsp_tag_self=4029;
p.bsp_tag_other=94;

p.bsp_params_file_name=['D:\Ayelet\Data\Data_Nlg_Proc\2bat_proj\bsp_params_bat_',num2str(bat),'.mat'];
p.path_day_dir=['D:\Ayelet\Data\2batproj\yr_2019_bat_',num2str(bat),'\',num2str(day)];
param_folder='D:\Ayelet\2bat_proj\Analysis\new_code\params\';
dir_param_file_name=fullfile(param_folder,'dirs_params.mat');
behav_param_file_name=fullfile(param_folder,'behav_params.mat');
load(dir_param_file_name)
p.sync_to='bsp';

p = BSP_PRE_PROC_data (p);
bsp_data =p.bsp_data;
tags = [p.bsp_data.tag_ID];
    tag_i = find(ismember(tags,p.bsp_tag_self));

behave_ts=[min(p.bsp_data.ts_us_upsampled) max(p.bsp_data.ts_us_upsampled)];
ball_pos_name=[ball_position_folder,'ball_pos_bat_',num2str(bat),'_day_',num2str(day),'.mat'];
[behavioral_modes]=find_behavioral_modes(bsp_data ,behav_param_file_name,tag_i,bat,day,behave_ts,dir_param_file_name,ball_pos_name); % assign bsp samples to different behaviors
