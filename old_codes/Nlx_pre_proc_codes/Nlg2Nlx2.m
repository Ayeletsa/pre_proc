
function Nlg2Nlx2(main_dir,post_rec_ref_channel)

% function Nlg2Nlx(main_dir)
%
%   main_dir: data dir
%
%
% Tamir Eliav
% revised - Didi Omer, May  2015
% revised - Didi Omer, June 2015
%

%This is Tamir�s function that Didi revised. This function create nlx files from the nlg files 
%both for the neural data and the events. 
%Note that this function will not create nlg files if there is already "nlg" folder for this 
%day in the relevant path. 
%% arrange files/folders

% case 1: no nlg and nlx directories


Nlg_InDir = fullfile(main_dir, 'nlg');
Nlx_OutDir = fullfile(main_dir, 'nlx');
if ~exist(Nlx_OutDir,'dir')
    mkdir(main_dir,'nlx')
else
    %     rmdir(Nlx_OutDir,'s');
    %     if ~exist(Nlx_OutDir,'dir')
    %         mkdir(main_dir,'nlx')
    %     end
    return
end
disp('1. constract/verify file structure...');
header_file = 'Nlg2Nlx_header.txt';

if ~exist(Nlg_InDir,'dir')
    mkdir(main_dir,'nlg');
    
    files = dir(main_dir);
    ii= find(~[files.isdir]==1);
    files_name = {};
    c = 1;
    for i=1:length(ii)
        [~,n,t] = fileparts(files(ii(i)).name);
        if strcmp(t,'.DAT') || ...
                strcmp(t,'.xlsx') || ...
                strcmp(t,'.CSV') || ...
                strcmp(t,'.NLE') ||...
                strcmp(t,'.txt'),
            
            files_name{c} = [n,t];
            c = c+1;
        end
    end
    
    
    
    if ~isempty(files_name)
        for i=1:length(ii)
            movefile([main_dir,'\',files(ii(i)).name],Nlg_InDir);
        end
    else
        if ~exist(Nlx_OutDir,'dir')
            warning('no Nlg files were found!');
            return;
        end
    end
    
else  % case the Nlg dir exists
    
    
    files1 = dir(main_dir);
    files2 = dir(Nlg_InDir);
    
    ii= find(~[files1.isdir]==1);
    files1_name = {};
    c = 1;
    for i=1:length(ii)
        [~,n,t] = fileparts(files1(ii(i)).name);
        if strcmp(t,'.DAT') || ...
                strcmp(t,'.xlsx') || ...
                strcmp(t,'.CSV') || ...
                strcmp(t,'.NLE') ||...
                strcmp(t,'.txt'),
            
            files1_name{c} = [n,t];
            c = c+1;
        end
    end
    
    ii= find(~[files2.isdir]==1);
    files2_name = {};
    c = 1;
    for i=1:length(ii)
        [~,n,t] = fileparts(files2(ii(i)).name);
        if strcmp(t,'.DAT') || ...
                strcmp(t,'.xlsx') || ...
                strcmp(t,'.CSV') || ...
                strcmp(t,'.NLE') ||...
                strcmp(t,'.txt'),
            
            files2_name{c} = [n,t];
            c = c+1;
        end
    end
    
    % case the Nlg dir exists but the files are outside
    if isempty(files2_name) && ~isempty(files1_name)
        for i=1:length(ii)
            movefile([main_dir,'\',files1(ii(i)).name],Nlg_InDir);
        end
    elseif isempty(files2_name) && isempty(files1_name)
        warning('no Nlg files were found!');
        return;
        
    end
    
end



%% parameters setting

motion_channel=1;
num_channels = 16;
data_cnl_ind = [0:15]; % note that this numbering system is of the neurologger which means channels 0-15
DATA_file_prefix = 'NEUR0';
% DATA_file_prefix = 'BACK';
zero_DC_level_bit = 2048;
is_invert_data = true;
is_remove_DC = true;
is_remove_flash_write_artifact = false;
use_clock_diff_correction = true;
% ref_channel 0 = use ground as refernce
if post_rec_ref_channel == 0
    use_post_rec_ref_channel = false;
else
    use_post_rec_ref_channel = true;
end

%% read EVENT file
event_file_name_xlsx = fullfile(Nlg_InDir, 'EVENTLOG.csv');
[NUM,TXT,RAW]=xlsread(event_file_name_xlsx);

% extract recording details from event file header
file_header_lines = TXT(1:3,1);
[splitstr] = regexp(file_header_lines{2}, '[;]+', 'split'); % 2nd header row
firmware_ver = regexp(splitstr{1}, '\d*[.]\d*','match');
serial_number = regexp(splitstr{2}, '\d*','match');
time = regexp(splitstr{3}, '\d*:\d*:\d*','match');
date = regexp(splitstr{4}, '\d*/\d*/\d*','match');
% ADC_period_usec = regexp(splitstr{5}, '\d*[.]\d*','match');
ADC_period_usec = regexp(splitstr{5}, '\d*','match');
[splitstr] = regexp(file_header_lines{3}, '[;]+', 'split'); % 3rd header row
ADC_resolution_uVolt = regexp(splitstr{1}, '\d*[.]\d*','match');

ADC_SAMPLE_PERIOD = str2num(cell2mat(ADC_period_usec))/num_channels*1e-6;
fs =  1/(ADC_SAMPLE_PERIOD * num_channels);
uVolt_per_bit = str2num(cell2mat(ADC_resolution_uVolt));
block_period_time_usec = (512/fs) * 1e6;
digital_data_block_period_time_usec = (512/fs) * 1e6;
block_period_time_usec_2 = 512 * (ADC_SAMPLE_PERIOD * num_channels) * 1e6;
file_len_time_usec = block_period_time_usec*1024; % TODO: calc
fs_acc = 1e6/(file_len_time_usec / 2048);
ts_offset_acc_neur_usec = file_len_time_usec / 2048;

%% extract events details
events_IX = NUM(1:end,1);
events_TS = NUM(1:end,3).*1e3;
events_TS_source = TXT(5:end,4);
events_type = TXT(5:end,5);
events_details = TXT(5:end,6);

%% join '...Continued' events
continued_event_lines_IX = find(isnan(events_IX))';
for ii_continued_event_line = continued_event_lines_IX
    % take the last line with a valid number in the event index column as
    % the event index
    last_valid_line = find( ~isnan(events_IX(1:ii_continued_event_line)), 1, 'last');
    events_details{last_valid_line} = [events_details{last_valid_line} events_details{ii_continued_event_line}];
end
events_IX(continued_event_lines_IX) = [];
events_TS(continued_event_lines_IX) = [];
events_TS_source(continued_event_lines_IX) = [];
events_type(continued_event_lines_IX) = [];
events_details(continued_event_lines_IX) = [];

%% join event type with event details to a single string
EventStrings = {};
for ii_event = 1:length(events_type)
    event_str = [events_type{ii_event} ' - ' events_details{ii_event}];
    EventStrings{ii_event} = event_str;
end

%% extract recording parameteres from the event file

disp('2. reading+saving parameters from event file...');

for i=1:size(events_type,1)
    if strcmp(events_type{i},'Recording parameters'),break;end % TODO: what is this line for?
end

paramst = strsplit(events_details{i},';');
for i=1:length(paramst)
    if length(paramst{i})>1
        t = strsplit(paramst{i},'=');
        t1 = t{1};
        t1= t1(find(~isspace(t1)));
        if strfind(t1,'-'), t1(find(t1=='-'))='';end
        if strfind(t1,'/'), t1(find(t1=='/'))='';end
        t{1} = t1;
        param.(t1) = t{2};
    end
end

% param.ChannelOverwrittenbyMotionSensor = str2double(char(regexp(param.ChannelOverwrittenbyMotionSensor,'(\d*)','match')));
% param.GyroscopeRangeIndex = str2double(char(regexp(param.GyroscopeRangeIndex,'(\d*)','match')));
% param.AccelerometerRangeIndex = str2double(char(regexp(param.AccelerometerRangeIndex,'(\d*)','match')));
param.LowCutoffFrequency = str2double(char(regexp(param.LowCutoffFrequency,'(\d*)','match')));

save( fullfile(Nlx_OutDir,'param.mat'), 'param' );


%% Apply clock difference correction (Transceiver vs. logger time)
% TODO: plot figure showing the clocks drift

% CD = clock difference
disp('3. correcting for clock diff...');

if use_clock_diff_correction
    
    PC_gen_events_IX = find(strcmp('PC-generated comment', events_type));
    CD_str = 'CD=';
    CD_values = [];
    CD_timestamps = [];
    for ii_event = 1:length(PC_gen_events_IX)
        curr_event_details = events_details{PC_gen_events_IX(ii_event)};
        CD_str_pos = strfind(curr_event_details, CD_str);
        if isempty(CD_str_pos)
            continue;
        end
        str = curr_event_details(CD_str_pos+3:end);
        CD = str2num(str(1:min(strfind(str,' '))));
        CD_values(end+1) = CD;
        CD_timestamps(end+1) = events_TS(PC_gen_events_IX(ii_event));
    end
    CD_values = CD_values.*1e6; % change from sec to usec
    
    % TODO: we need to check what is the meaning of CD (clock difference). Logger-Tx
    % or Tx-Logger?
    % Is this old comment relevant? Didi
    % TODO: best to do is to check this again after Jacob bring us the new
    % version that fix the problem from version 1.69
    
    CD_event_ts__logger_time = CD_timestamps-CD_timestamps(1);        % removed 1st timestamp to balance the fit data closer to zero (otheriwse it gives a warning message)
    CD_event_ts__Tx_time = CD_event_ts__logger_time - CD_values;
    transceiver_2_logger_time_fit = polyfit(CD_event_ts__Tx_time, CD_event_ts__logger_time , 1);
    ts_source_transceiver_events_IX = find(strcmp('Transceiver', events_TS_source));
    events_TS(ts_source_transceiver_events_IX) = polyval(transceiver_2_logger_time_fit , events_TS(ts_source_transceiver_events_IX) - CD_timestamps(1)) + CD_timestamps(1);
    for ii_event = 1:length(ts_source_transceiver_events_IX)
        events_TS_source{ts_source_transceiver_events_IX(ii_event)} = 'Logger';
    end
    
    % TODO: plot some figure, and make a check that the clock dcorrection
    % was OK
end

%% write event file in Nlx format

disp('4. writing event file...');
NlxEventFile = fullfile(Nlx_OutDir, 'EVENTS.nev');
Mat2NlxEV(NlxEventFile, 1, 1, [], [1 0 0 0 1 0] , events_TS', EventStrings');

%% create sperate event file for each event category
event_type_list = unique(events_type);
for ii_event_type = 1:length( event_type_list )
    event_type_string = event_type_list{ii_event_type};
    event_type_events_IX = find(strcmp( event_type_string , events_type));
    Nlx_event_type_file = fullfile(Nlx_OutDir, ['EVENTS__' event_type_string '.nev']);
    Mat2NlxEV(Nlx_event_type_file , 1, 1, [], [1 0 0 0 1 0] , events_TS(event_type_events_IX)', EventStrings(event_type_events_IX)');
end

%% create battery discharge plot
disp('5. calculating battery discharge...');
figure;
PC_gen_events_IX = find(strcmp('PC-generated comment', events_type));
mode_changed_events_IX = find(strcmp('Mode change', events_type));
BV_str = 'BV=';
BV_values = [];
BV_timestamps = [];
for ii_event = 1:length(PC_gen_events_IX)
    curr_event_details = events_details{PC_gen_events_IX(ii_event)};
    BV_str_pos = strfind(curr_event_details, BV_str);
    if isempty(BV_str_pos)
        continue;
    end
    
    BV_values_str_interval = (BV_str_pos:BV_str_pos+4) + length(BV_str);
    BV_values(end+1) = str2num( curr_event_details(BV_values_str_interval) );
    BV_timestamps(end+1) = events_TS(PC_gen_events_IX(ii_event));
end
usec_2_min = 1/60*1e-6;
% BV_timestamps = BV_timestamps(4:end);
% BV_values = BV_values(4:end)
plot( usec_2_min.*(BV_timestamps - BV_timestamps(1)), BV_values , '-o')
title('Battery discharge')
xlabel('time (minutes)')
ylabel('Voltage (V)')
grid on

hold on
mode_changed_events_ts = events_TS(mode_changed_events_IX) - BV_timestamps(1);
ylimits = get(gca, 'ylim');
plot( repmat(usec_2_min .* mode_changed_events_ts,1,2)' , ylimits', '--m' )

legend({'battey voltage','record mode change event'})

fig_file = fullfile(Nlx_OutDir, 'Battery_discharge');
saveas(gcf, fig_file , 'jpg')
saveas(gcf, fig_file , 'fig')
close all

%% create Nlx data files - continuous sampling files called (.ncs)

% identify  'File started' events and take timestamps
FileStarted_IX = strcmp('File started', events_type);
FileStarted_TS = events_TS(FileStarted_IX);
FileStarted_details = events_details(FileStarted_IX);

% read header template
header = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');
% change fields that are recording specific
header_neural = strrep(header,'-SamplingFrequency', ['-SamplingFrequency ' char(vpa(fs))]);
header_motion = strrep(header,'-SamplingFrequency', ['-SamplingFrequency ' char(vpa(fs_acc))]);
disp('6. Nlg -> Nlx...');

for ii_file_start_entry = 1:length(FileStarted_TS)
    
    %%
    file_str = FileStarted_details{ii_file_start_entry};
    file_name = [DATA_file_prefix , file_str(end-2:end) '.DAT'];
    %     file_name = [DATA_file_prefix '_' file_str(end-2:end) '.DAT'];
    fid = fopen(fullfile(Nlg_InDir, file_name));
    filedata = fread(fid, 'uint16', 0, 'l');
    fclose(fid);
    filedata = double(filedata);
    data = reshape(filedata, num_channels, length(filedata)/num_channels);
    file_TS_usec = FileStarted_TS(ii_file_start_entry);
    
    
    if isfield(param, 'MotionSensorLogging');
        if strcmp(param.MotionSensorLogging,' Enabled')
            
            
            %% parse accelerometer data (digital)
            block_data = data(param.ChannelOverwrittenbyMotionSensor+1,:);
            block_data = reshape(block_data,256,2048);
            acc_data = [];
            gyro_data = [];
            magnet_data = [];
            digital_data_ts = [];
            for ii_rec = 1:2048
                %% check barker
                rec_data = block_data(:,ii_rec);
                barker1 = rec_data(1);
                barker2 = rec_data(2);
                if (barker1~=13579) || (barker2~=24680)
                    acc_data    = [acc_data,    [0 0 0]'];
                    gyro_data   = [gyro_data,   [0 0 0]'];
                    magnet_data = [magnet_data, [0 0 0]'];
                    warning('wrong BARKER')
                    continue;
                end
                acc_offset = rec_data(3);
                gyro_offset = rec_data(4);
                magnet_offset = rec_data(5);
                acc_len = rec_data(7);
                gyro_len = rec_data(8);
                magnet_len = rec_data(9);
                
                
                % Didi - for each of the 2048 sequences take the mean for each sensor.
                
                acc_data_seg = rec_data( (acc_offset+1) : (acc_offset+acc_len));
                gyro_data_seg = rec_data( (gyro_offset+1) : (gyro_offset+gyro_len));
                magnet_data_seg = rec_data( (magnet_offset+1) : (magnet_offset+magnet_len));
                if ~isempty(acc_data_seg)
                    acc_data_seg = acc_data_seg(1:floor(length(acc_data_seg)/3)*3);
                    gyro_data_seg = gyro_data_seg(1:floor(length(gyro_data_seg)/3)*3);
                    magnet_data_seg = magnet_data_seg(1:floor(length(magnet_data_seg)/3)*3);
                    
                    acc_data_seg = median(reshape(acc_data_seg,[3,length(acc_data_seg)/3]),2);
                    gyro_data_seg = median(reshape(gyro_data_seg,[3,length(gyro_data_seg)/3]),2);
                    magnet_data_seg = median(reshape(magnet_data_seg,[3,length(magnet_data_seg)/3]),2);
                    
                    acc_data    = [acc_data,    acc_data_seg];
                    gyro_data   = [gyro_data,   gyro_data_seg];
                    magnet_data = [magnet_data, magnet_data_seg];
                    
                else
                    % Didi - case the motion data is empty (first sample, first
                    % trial
                    acc_data    = [acc_data,    [0;0;0]];
                    gyro_data   = [gyro_data,   [0;0;0]];
                    magnet_data = [magnet_data, [0;0;0]];
                    
                end
                
                
                
                
            end
            %     digital_data_ts = digital_data_ts';
            %     digital_data_ts_usec = digital_data_ts .* 1e3;
            
            
            
            
            
            
            acc_data_X = acc_data(1,:);
            acc_data_Y = acc_data(2,:);
            acc_data_Z = acc_data(3,:);
            gyro_data_X = gyro_data(1,:);
            gyro_data_Y = gyro_data(2,:);
            gyro_data_Z = gyro_data(3,:);
            magnet_data_X = magnet_data(1,:);
            magnet_data_Y = magnet_data(2,:);
            magnet_data_Z = magnet_data(3,:);
            
            % convert to signed values
            uint16_to_int16 = @(a)(a.*(a<=(2^15-1)) + (a-2^16).*(a>(2^15-1)));
            acc_data_X = uint16_to_int16(acc_data_X);
            acc_data_Y = uint16_to_int16(acc_data_Y);
            acc_data_Z = uint16_to_int16(acc_data_Z);
            gyro_data_X = uint16_to_int16(gyro_data_X);
            gyro_data_Y = uint16_to_int16(gyro_data_Y);
            gyro_data_Z = uint16_to_int16(gyro_data_Z);
            magnet_data_X = uint16_to_int16(magnet_data_X);
            magnet_data_Y = uint16_to_int16(magnet_data_Y);
            magnet_data_Z = uint16_to_int16(magnet_data_Z);
            
            acc_data_X_blocks = reshape(acc_data_X,[],4)';
            acc_data_Y_blocks = reshape(acc_data_Y,[],4)';
            acc_data_Z_blocks = reshape(acc_data_Z,[],4)';
            gyro_data_X_blocks = reshape(gyro_data_X,[],4)';
            gyro_data_Y_blocks = reshape(gyro_data_Y,[],4)';
            gyro_data_Z_blocks = reshape(gyro_data_Z,[],4)';
            magnet_data_X_blocks = reshape(magnet_data_X,[],4)';
            magnet_data_Y_blocks = reshape(magnet_data_Y,[],4)';
            magnet_data_Z_blocks = reshape(magnet_data_Z,[],4)';
            
            temp = linspace(0,file_len_time_usec,2049);
            temp = temp(1:end-1);
            digital_data_ts_usec = file_TS_usec + temp;
            digital_data_ts_usec = digital_data_ts_usec - ts_offset_acc_neur_usec;
            digital_data_block_ts_usec = digital_data_ts_usec(1:512:end);
            digital_data_block_ts_usec =digital_data_block_ts_usec+(8.7381333/1000);
            digital_data_block_ts_usec = round(digital_data_block_ts_usec);          % usec is good enough... (and we can have double/integer problems...)
            
            acc_parsed_data = cat(3,...
                acc_data_X_blocks,...
                acc_data_Y_blocks,...
                acc_data_Z_blocks,...
                gyro_data_X_blocks,...
                gyro_data_Y_blocks,...
                gyro_data_Z_blocks,...
                magnet_data_X_blocks,...
                magnet_data_Y_blocks,...
                magnet_data_Z_blocks...
                );
            
            acc_ch_names = {...
                'acc_X',...
                'acc_Y',...
                'acc_Z',...
                'gyro_X',...
                'gyro_Y',...
                'gyro_Z',...
                'magnet_X',...
                'magnet_Y',...
                'magnet_Z',...
                };
            
            % write motion data to nlx file format
            for ii_acc_chnl = 1:9
                file_name = ['Motion_' acc_ch_names{ii_acc_chnl} '.ncs'];
                out_file = fullfile(Nlx_OutDir, file_name);
                
                % we do this because Nlx functions have bugs with writing header
                % when using append... (this way works...)
                if exist(out_file, 'file')
                    append_flag = 1;
                else
                    append_flag = 0;
                end
                
                blocks_to_write = squeeze(acc_parsed_data(:,:,ii_acc_chnl))';
                Mat2NlxCSC(out_file, append_flag, 1, 0, [1 0 1 0 1 1], digital_data_block_ts_usec, repmat(fs_acc,1,4), blocks_to_write , header_motion );
                
                %             for ii_block = 2:2
                %                 data_chunk = squeeze(acc_parsed_data(ii_block,:,ii_acc_chnl))';
                %                 Mat2NlxCSC(out_file, append_flag, 1, 0, [1 0 0 0 1 1], digital_data_block_ts_usec(ii_block), data_chunk, header_motion );
                %             end
            end
            
        end
    end
    
    %% remove DC
    if is_remove_DC
        for ii_channel = 1:size(data,1)
            %         data(ii_channel,:) = data(ii_channel,:) - mean(data(ii_channel,:));  % TAMIR: we prefer not to use this since data could be biased around real zero voltage (artifacts for example...)
            data(ii_channel,:) = data(ii_channel,:) - zero_DC_level_bit;         % There should be a theoretical constant representing the real zero voltage, but because the INTAN have a DC dependency on freq it is not perfect
        end
    end
    
    %% Remove repetative flash write artifact
    % % % %     if is_remove_flash_write_artifact
    % % % %         weak_artifacts_block_IX = [1:256 513:768 1025:1280 1537:1793];
    % % % %         strong_artifacts_block_IX = setdiff(1:2048, weak_artifacts_block_IX);
    % % % %         trim_prc = 1;
    % % % %         for ii_channel = 1:size(data,1)
    % % % %             data_channel_blocks = reshape( data(ii_channel,:), 256, [] );
    % % % %             flash_write_artifact_shape_weak = trimmean(data_channel_blocks(:,weak_artifacts_block_IX), trim_prc , 'round', 2);
    % % % %             flash_write_artifact_shape_strong = trimmean(data_channel_blocks(:,strong_artifacts_block_IX), trim_prc , 'round', 2);
    % % % %
    % % % %             %% TEMP
    % % % %             sdf = data_channel_blocks(:,strong_artifacts_block_IX);
    % % % %             sdf = sdf - mean(mean(sdf));
    % % % %             sdf = abs(sdf);
    % % % %             sdf = sdf.*uVolt_per_bit;
    % % % %             figure
    % % % %             hold all
    % % % %             plot(mean(sdf,2));
    % % % %             plot(prctile(sdf,50,2));
    % % % %             plot(prctile(sdf,60,2));
    % % % %             plot(prctile(sdf,70,2));
    % % % %             plot(prctile(sdf,80,2));
    % % % %             plot(prctile(sdf,90,2));
    % % % %             legend({'averaged','50%','60%','70%','80%','90%'});
    % % % %             ylim([0 50])
    % % % %
    % % % %             %%
    % % % %             data_channel_blocks_corrected = zeros(size(data_channel_blocks));
    % % % %             for ii_block = weak_artifacts_block_IX
    % % % %                 data_channel_blocks_corrected(:,ii_block) = data_channel_blocks(:,ii_block) - flash_write_artifact_shape_weak;
    % % % %             end
    % % % %             for ii_block = strong_artifacts_block_IX
    % % % %                 data_channel_blocks_corrected(:,ii_block) = data_channel_blocks(:,ii_block) - flash_write_artifact_shape_strong;
    % % % %             end
    % % % %             data(ii_channel,:) = reshape(data_channel_blocks_corrected,1,[]);
    % % % %        end
    % % % %     end
    
    %% if we recorded with GND as ref. channel we want to have a
    % "post-recording ref. channel"
    if use_post_rec_ref_channel
        gain_relative_to_ref_chnl = 1.*ones(1,16);
        for ii_channel = 1:size(data,1)
            if ii_channel == post_rec_ref_channel
                continue;
            end
            data(ii_channel,:) = data(ii_channel,:) - gain_relative_to_ref_chnl(ii_channel).*data(post_rec_ref_channel,:);
        end
    end
    
    %% change to uVolt units and invert
    data = data.*uVolt_per_bit;
    if is_invert_data
        % invert data
        data = -data;
    end
    
    
    %% write data chunks to Nlx file format
    %     file_TS_usec = FileStarted_TS(ii_file_start_entry);
    disp(['writing file ' file_str])
    for cnl = data_cnl_ind
        cnl_data = data(cnl+1,:);
        cnl_data_blocks = vec2mat(cnl_data,512)';
        num_blocks = size(cnl_data_blocks, 2);
        blocks_timestamps_usec = file_TS_usec + (0:block_period_time_usec:(block_period_time_usec*(num_blocks-1)));
        file_name = ['CSC' num2str(cnl) '.ncs'];
        out_file = fullfile(Nlx_OutDir, file_name);
        
        %         we do this because Nlx functions have bugs with writing header
        %         when using append... (this way works...)
        if exist(out_file, 'file')
            append_flag = 1;
        else
            append_flag = 0;
        end
        %append_flag
        
        Mat2NlxCSC(out_file, append_flag, 1, 0, [1 0 0 0 1 1], blocks_timestamps_usec, cnl_data_blocks, header_neural );
    end
    %pause;
end

% %
% % %% try to add header
% % file = 'D:\arseny\NlgRec\neurologger_recording_20140107_arseny\test\test.ncs';
% % % Mat2NlxCSC(file, 0, 1, 0, [1 0 0 0 1 1], zeros(1,0), zeros(512,0), header );
% % % Mat2NlxCSC(file, 0, 1, 0, [1 0 0 0 1 1], blocks_timestamps_usec, cnl_data_blocks, header );
% % Mat2NlxCSC(file, 1, 1, 0, [1 0 0 0 1 1], blocks_timestamps_usec+20*1e6, cnl_data_blocks, header );

%%
disp('Finished converting all files from Nlg format to Nlx format!')









%%