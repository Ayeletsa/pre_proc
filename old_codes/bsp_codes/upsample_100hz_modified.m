function bsp_data_new = upsample_100hz_modified(bsp_data)
%bsp_data - bsp struct after linearization of the data, in ~18 Hz sample
%rate. 
%bsp_data - new struct of 100Hz saple rate of the sampe data. including a
%new field: non_spikes_ts regarding gaps in which interpolation was not
%performed due to a bspoon problem (no samples in over ~2 sec). the spikes
%in these times will not be onsidered in place fields analysis. 
bsp_data_new = bsp_data;
bsp_freq = 1/(min(diff(bsp_data(1).ts_ns))*10^-9); T = (1/bsp_freq)*10^9; %musec
bsp_new_freq = 100; T_new = (1/bsp_new_freq)*10^9; %musec
min_sec_to_interp=2;

%find time vector for both
min_time=max([min(bsp_data(1).ts_ns), min(bsp_data(2).ts_ns)]);
max_time=min([max(bsp_data(1).ts_ns), max(bsp_data(2).ts_ns)]);

%new time vector

 for tag_i=1:length(bsp_data)  
     
t_new_vec=min_time:T_new:max_time;
new_pos=NaN*zeros(length(t_new_vec),3);

ts_ns = bsp_data(tag_i).ts_ns;
% if min(ts_ns)<min_time 
%     ts_ns(find(ts_ns<min_time))=[];
%     
% end
% 
% if max(ts_ns)<max_time
%     ts_ns(find(ts_ns>max_time))=[];
% end

diff_ts_ns = diff(ts_ns);
bad_ind = find(diff_ts_ns>min_sec_to_interp*bsp_freq*T); %hold the samples that have a gap to the 
% near sample of more than ~2sec.
blocks_to_interpolate=[];
non_spikes_ts=[];
if isempty(bad_ind)
    blocks_to_interpolate = [1 length(ts_ns)];
    non_spikes_ts = [];
else
    blocks_to_interpolate(1,:) = [1 bad_ind(1)];
    non_spikes_ts = [];
    for i = 1:length(bad_ind)-1
        if bad_ind(i)+1 == bad_ind(i+1)
            continue
        end
        blocks_to_interpolate(end+1,:) = [(bad_ind(i)+1)  bad_ind(i+1)];
        non_spikes_ts(end+1,:) = [ts_ns(blocks_to_interpolate(end-1,2)) ts_ns(blocks_to_interpolate(end,1))];
    end
    blocks_to_interpolate(end+1,:) = [(bad_ind(end)+1)  length(ts_ns)];
    non_spikes_ts(end+1,:) = [ts_ns(blocks_to_interpolate(end-1,2)) ts_ns(blocks_to_interpolate(end,1))];
end
    
 pos = bsp_data(tag_i).normalized_pos;
for i=1:length(blocks_to_interpolate(:,1))
    %find the closest index in the new time vector:
    [~, closest_start]=min(abs(ts_ns(blocks_to_interpolate(i,1))-t_new_vec));
     [~, closest_end]=min(abs(ts_ns(blocks_to_interpolate(i,2))-t_new_vec));

    
    %use the new time vector for interpolation:    
    new_ts = t_new_vec(closest_start:closest_end);

   % interpolate and use it as the new position: 
    new_pos(closest_start:closest_end,:) = interp1(ts_ns,pos,new_ts','linear');
end

% pos = bsp_data(tag_i).pos;
% new_pos = interp1(ts_ns,pos,new_ts,'linear');


bsp_data_new(tag_i).pos_linear_upsamp= new_pos;
bsp_data_new(tag_i).ts_usec_upsamp = t_new_vec';
bsp_data_new(tag_i).non_spikes_ts = non_spikes_ts; %time stamps that should not be considered for spikes.

 end