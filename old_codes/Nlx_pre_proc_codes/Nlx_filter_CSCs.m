function  Nlx_filter_CSCs(expname,forcecalc)

%
% Didi Omer September 2015 

%%
eval(expname); 

%%
fwin = 2;         % we will run over the data in 2-min windows but save           % Same as we use for filtering ripples
passband_spikes   = [600 6000];        % Filter for spikes
passband_LFP      = [0.5 400];         % Filter for LFPs
Input_Dir_Recording = param.path.Nlx;
times_microsec_total = param.spikes.neural_data_ts;
csc_dir_LFP    = param.path.CSC_LFP;
csc_dir_spikes = param.path.CSC_spikes;

%% extract LFPs and save them
% skip CSC extraction and exit function if we have already extracted the data
run_LFP_filtering = 0;
if ~exist(csc_dir_LFP,'dir')
    mkdir(csc_dir_LFP);
    run_LFP_filtering = 1;
end
if forcecalc
    % TODO: delete/rename old files!!!
    run_LFP_filtering = 1;
end
if ~param.spikes.generate_LFP_files
    run_LFP_filtering = 0;
end
if run_LFP_filtering
    
    fprintf('Extracting LFP from Tetrode(16): 0000');
    for ii_channel=1:(length(param.tetrodes.TT)*4)

        filename_in = fullfile(Input_Dir_Recording,['CSC',num2str(ii_channel-1),'.ncs']);
        filename_out = fullfile(csc_dir_LFP, ['LFP_', param.expname, '_TT', num2str(ceil(ii_channel/4)), '_ch' num2str(mod(ii_channel-1,4))]);
        [signal,ts,fs]= FIRfilterNlxSlide(filename_in, times_microsec_total, passband_LFP ,2);

        % decimate LFP
%         subsample = ceil(1000/passband_LFP(2)/(1000/fs) / 2);
        subsample = 10;
        fs_LFP = fs / subsample ;
%         fs_LFP = 1000/(subsample*(1000/fs));
        signal = signal(1:subsample:end);
        ts = ts(1:subsample:end);

        LFP.params.fwin =  fwin;
        LFP.params.passband = passband_LFP ;
        LFP.params.CSC_SamplingFreq = fs_LFP;
        LFP.params.times_microsec_total = times_microsec_total;
        LFP.params.Recording_directory = Input_Dir_Recording;
    %     LFP.data.Samples = signal ;
    %     LFP.data.Timestamps = ts ;
        save([filename_out '.mat'],'LFP');
        nlx_csc_write([filename_out '.ncs'], signal, ts, fs_LFP);
        clear LFP sig x
        fprintf('\b\b\b\b%04d%', ii_channel);
    end
    fprintf('\n');

else
    disp('======================');
    disp (['Skipping LFP extraction. Already extracted in ' csc_dir_LFP]);
    disp('======================');
end

%% exctract high-pass for spike detection

run_SPIKES_filtering = 0;
% skip CSC extraction and exit function if we have already extracted the data
if isempty (dir(csc_dir_spikes))
    mkdir(csc_dir_spikes);
    run_SPIKES_filtering = 1;
end
if forcecalc
    % TODO: delete/rename old files!!!
    run_SPIKES_filtering = 1;
end

if run_SPIKES_filtering
    
    fprintf('Extracting high-pass CSC from Tetrode(16): 0000');
    for ii_channel=1:(length(param.tetrodes.TT)*4)

        filename_in = fullfile(Input_Dir_Recording,['CSC',num2str(ii_channel-1),'.ncs']);
        filename_out = fullfile(csc_dir_spikes, ['spikes_', param.expname, '_TT', num2str(ceil(ii_channel/4)), '_ch' num2str(mod(ii_channel-1,4))]);

        [signal,ts,fs]= FIRfilterNlxSlide(filename_in, times_microsec_total, passband_spikes, 2);

        Spikes.params.fwin =  fwin;
        Spikes.params.passband = passband_spikes ;
        Spikes.params.CSC_SamplingFreq = fs;
        Spikes.params.times_microsec_total = times_microsec_total;
        Spikes.params.Recording_directory = Input_Dir_Recording;
    %     Spikes.data.Samples = signal;
    %     Spikes.data.Timestamps = ts;
        save(filename_out,'Spikes');
        nlx_csc_write([filename_out '.ncs'], signal, ts, fs);
        clear signal ts fs
        fprintf('\b\b\b\b%04d%', ii_channel);
    end
    fprintf('\n');
    
else
    disp('======================');
    disp (['Skipping high-pass extraction. Already extracted in ' csc_dir_spikes]);
    disp('======================');
    return;
end;



