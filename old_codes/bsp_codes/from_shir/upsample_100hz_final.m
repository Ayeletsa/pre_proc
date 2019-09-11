function bsp_pos_new = upsample_100hz_final(bsp_pos)
%bsp_pos - bsp struct after linearization of the data, and after hole
%filling. in holes which should not be intewrepolated, there should be nans
%in the position vector.
%bsp_pos_new struct of 100Hz saple rate of the sampe data. including a
%
%new field: non_spikes_ts regarding gaps in which interpolation was not
%performed due to a bspoon problem. the spikes in these times will not be 
%considered in place fields analysis. 

bsp_new_freq = 100; T_new = (1/bsp_new_freq)*10^6; %musec

ts_nlg = bsp_pos.ts_nlg_usec;
pos = bsp_pos.pos_linear;
new_ts = ts_nlg(1) : T_new : ts_nlg(end);

new_pos = interp1(ts_nlg,pos,new_ts,'linear');

nan_pos = isnan(new_pos);
tmp = diff([0 nan_pos]);
bad_ind_start = find(tmp==1)-1;
bad_ind_end = find(tmp==-1);
non_spikes_ts = [new_ts(bad_ind_start)' new_ts(bad_ind_end)'];

bsp_pos_new = bsp_pos;
bsp_pos_new.pos_linear = new_pos';
bsp_pos_new.ts_nlg_usec = new_ts';
bsp_pos_new.non_spikes_ts = non_spikes_ts; %time stamps that should not be considered for spikes.





