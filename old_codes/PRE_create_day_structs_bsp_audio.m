function PRE_create_day_structs_bsp_audio
folder_name_bat='D:\Ayelet\Data\2batproj\yr_2018_bat_2389\20180815\aud_2389\';
folder_name_other='D:\Ayelet\Data\2batproj\yr_2018_bat_2389\20180815\aud_2287\';

p.path_day_dir='D:\Ayelet\Data\2batproj\yr_2018_bat_2389\20180815';
p.bsp_tag_self=484;
p.bsp_tag_other=4029;

%% 1. nlg 2 nlx

nlg2nlx_audio(folder_name_bat)
nlg2nlx_audio(folder_name_other)

%% 2. create audio file
Filename=[folder_name_bat,'audio.ncs'];
[signal, ts, fs] = Nlx_csc_read(Filename,[]);
p.audio.siganl_bat=signal;
p.audio.ts_bat=ts;
p.audio.fs_bat=fs;

Filename=[folder_name_other,'audio.ncs'];
[signal, ts, fs] = Nlx_csc_read(Filename,[]);
p.audio.siganl_other=signal;
p.audio.ts_other=ts;
p.audio.fs_other=fs;

%% 3. extract BSP
p = BSP_PRE_extract_data (p); %TO DO - add landmarks  + take parameters to the param_in script

%% 4. sync
% need to work on it...
p= PRE_sync_nlg2bsp(p);
p= PRE_sync_nlg2bsp(p);


   
%a. for self
ts_bat_msec   = p.audio.ts_bat   / 1e3;
p.ts_bat_msec_fitted_to_bsp_msec = polyval (p.time_conv_p_msec_nlx2bsp_self,ts_bat_msec);

%b. for other
ts_other_msec   = p.audio.ts_other   / 1e3;
p.ts_other_msec_fitted_to_bsp_msec = polyval (p.time_conv_p_msec_nlx2bsp,ts_other_msec);


%% 5. save

struct_name='D:\Ayelet\Data\Data_Nlg_Proc\yr_2018_bat_2389\aud_day_structs\day_20180815.mat';
save(struct_name,'p','-v7.3')

 


end
