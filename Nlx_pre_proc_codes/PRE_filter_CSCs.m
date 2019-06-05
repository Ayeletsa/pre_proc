function PRE_filter_CSCs(p)

main_dir=p.path_day_dir;
active_channels=p.active_channels;
bat=p.bat;
path_dataout=p.path_dataout;
year_bat_path=p.year_bat_path;

main_dir_in = main_dir;
date = char (regexp (main_dir,'\d{8}','match'));
bat_id = num2str(bat);
main_dir_out = [path_dataout,'\',year_bat_path,'\',date];

if exist (main_dir_out,'dir')
    return
end

%% params
load(p.filter_params_file_name);

active_channels = reshape(active_channels',[],1)';

%%
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end
if poolsize ~=6  % checking to see if my pool is already open
    parpool(6)
else
    disp('matlab workers already open')
end

%% extract LFPs and save them
% skip CSC extraction and exit function if we have already extracted the data


t_start_end = [];
clear filter_params
filter_params.passband  = passband_LFP;
filter_params.fwin      = fwin;
filter_params.resample_fs = LFP_resamlpe_fs;

LFP_out_dir = fullfile(main_dir_out,'LFP');
if ~exist(LFP_out_dir,'dir')
    mkdir(LFP_out_dir)
    
    parfor ii_ch = 1:length(active_channels)
        if ~active_channels(ii_ch)
            continue;
        end
        TT = ceil(ii_ch/4);
        ch_num = mod(ii_ch-1,4)+1;
        file_IN = fullfile(main_dir_in,'nlx',['CSC' num2str(ii_ch-1) '.ncs']);
        file_OUT = fullfile(LFP_out_dir,['LFP_bat_',bat_id,'_day_',date,'_TT_' num2str(TT) '_ch' num2str(ch_num) '.ncs'])
        
        Nlx_filter_CSC2_tamir(file_IN, file_OUT, t_start_end, filter_params)
    end
    
else
    disp('======================');
    disp (['Skipping LFP extraction. Already extracted in ' fullfile(main_dir_out,'LFP')]);
    disp('======================');
end

%% exctract high-pass for spike detection

t_start_end = [];
clear filter_params
filter_params.type = spikes_filter_type;
filter_params.passband  = passband_spikes(1);

filter_params.fwin      = fwin;

spikes_out_dir = fullfile(main_dir_out,'spikes');
if ~exist(spikes_out_dir,'dir')
    mkdir(spikes_out_dir)
    
    parfor ii_ch = 1:length(active_channels)
        if ~active_channels(ii_ch)
            continue;
        end
        TT = ceil(ii_ch/4);
        ch_num = mod(ii_ch-1,4)+1;
        file_IN = fullfile(main_dir_in,'nlx',['CSC' num2str(ii_ch-1) '.ncs']);
        file_OUT = fullfile(spikes_out_dir,['spikes_bat_',bat_id,'_day_',date,'_TT' num2str(TT) '_ch' num2str(ch_num) '.ncs'])
        
        Nlx_filter_CSC2_tamir(file_IN, file_OUT, t_start_end, filter_params)
    end
    
else
    disp('======================');
    disp (['Skipping high-pass extraction. Already extracted in ' fullfile(main_dir_out,'spikes')]);
    disp('======================');
    return;
end


%%

poolobj = gcp('nocreate');
delete(poolobj);


end