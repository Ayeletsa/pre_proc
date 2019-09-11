function new_linearize_bsp_pos = linearize_bsp(bsp_pos_struct,tunnel_calib_path)
% 
% bsp_pos_struct : a struct with the raw bsp data that have the following field: bsp_pos.pos
% tunnel_calib_path is the path for the relevant tunnel calibration - maybe
% add this in future to the param of each experiment. for example:
%'D:\Shir\Project\Tamir\yr2017_bat0148_behav_neural\calibrations\20170809__tunnel_midline\calib_tunnel.mat'

dis_TH = 2; %m this is the maximum distance of a single sample, away from 
% the midline of the tunnel that we allow.

load(tunnel_calib_path);
[bsp_pos_linearized,dis,t] = distance2curve(calib_tunnel.curvexy,bsp_pos_struct.pos(:,1:2),'linear');

%Remove outlier samples:
bad_ind = find(dis>=dis_TH);
bsp_pos_linearized(bad_ind) = [];
t(bad_ind) = [];
bsp_pos_struct.ts_nlg_usec(bad_ind) = [];

bsp_pos_struct.pos_linear = t*calib_tunnel.tunnel_length;

new_linearize_bsp_pos = bsp_pos_struct;

  