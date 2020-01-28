function Nlg2Nlx_64ch(main_dir)

% function Nlg2Nlx(main_dir)
%
%   main_dir: data dir
%
%
% Tamir Eliav
% revised - Didi Omer, May  2015
% revised - Didi Omer, June 2015
%


%% arrange files/folders

disp('1. constract/verify file structure...');
header_file = 'Nlg2Nlx_header.txt';

% case 1: no nlg and nlx directories


Nlg_InDir = fullfile(main_dir, 'nlg');
Nlx_OutDir = fullfile(main_dir, 'nlx');
if ~exist(Nlg_InDir,'dir')
    error('NLG input folder does not exist');
end

if exist(Nlx_OutDir,'dir')
    return
end
if ~exist(Nlx_OutDir,'dir')
    mkdir(main_dir,'nlx')
end
%% parameters setting
num_channels = 64;
data_cnl_ind = [0:63]; % note that this numbering system is of the neurologger which means channels 0-15
DATA_file_prefix = 'NEUR';
zero_DC_level_bit = 2^15;
is_invert_data = true;
is_remove_DC = true;
is_remove_flash_write_artifact = false;
use_clock_diff_correction = false;
use_post_rec_ref_channel = false; % use this if you recorded with GND as ref channel
post_rec_ref_channel = 1; % (1-16) If the above is true - choose the channel you want to substract from all the other channels


%% Neurologger mapping to EIB mapping:
% assuming the logger is connected to the EIB so the SD card is in the side
% of TT1 and TT16 in the EIB.
ch_per_TT_name = zeros([1,num_channels]);
TT_name = [1 9 1 9 1 9 1 9 2 10 2 10 2 10 2 10 7 15 7 15 7 15 7 15 ...
    8 16 8 16 8 16 8 16 3 11 3 11 3 11 3 11 4 12 4 12 4 12 4 12 ...
    5 13 5 13 5 13 5 13 6 14 6 14 6 14 6 14];
tmp = (ones([num_channels/8, 1])*[1:4])';
ch_per_TT_name(1:2:end) = tmp(:);
ch_per_TT_name(2:2:end) = tmp(:);

%% read EVENT file
event_file_name_xlsx = fullfile(Nlg_InDir, 'EVENTLOG.csv');
[NUM,TXT,RAW]=xlsread(event_file_name_xlsx);

% extract recording details from event file header
file_header_lines = TXT(1,1);
[splitstr] = regexp(file_header_lines{1}, '[;]+', 'split'); % 2nd header row
firmware_ver = regexp(splitstr{1}, '\d*[.]\d*','match');
serial_number = regexp(splitstr{2}, '\d*','match');
time = regexp(splitstr{3}, '\d*:\d*:\d*','match');
date = regexp(splitstr{4}, '\d*/\d*/\d*','match');
ADC_period_usec = regexp(splitstr{5}, '\d*','match');
ADC_resolution_uVolt = regexp(splitstr{6}, '\d*[.]\d*','match');

samples_per_block = 512; %64-ch logger
blocks_per_file = 256; %64-ch logger
ADC_SAMPLE_PERIOD = str2num(cell2mat(ADC_period_usec))/num_channels*1e-6;%sec
fs = 32000; %1/(ADC_SAMPLE_PERIOD * num_channels);%Hz
uVolt_per_bit = str2num(cell2mat(ADC_resolution_uVolt));
block_period_time_usec = (512/fs) * 1e6;
digital_data_block_period_time_usec = (512/fs) * 1e6;
block_period_time_usec_2 = 512 * (ADC_SAMPLE_PERIOD * num_channels) * 1e6;
% file_len_time_usec = block_period_time_usec*1024; % TODO: calc
file_len_time_usec = block_period_time_usec*256;
fs_acc = 1e6/(file_len_time_usec / 2048);
ts_offset_acc_neur_usec = file_len_time_usec / 2048;

%% extract events details
events_IX = NUM(1:end,1);
events_TS = NUM(1:end,3).*1e3;
events_TS_source = TXT(3:end,4);
events_type = TXT(3:end,5);
events_details = TXT(3:end,6);

%% join '...Continued' events
continued_event_lines_IX = find(strcmp(events_type,'...Continued'));
for ii_continued_event = 1:length(continued_event_lines_IX)
    % take the last line with a valid number in the event index column as
    % the event index
    last_valid_line = find( ~strcmp(events_type(1:continued_event_lines_IX(ii_continued_event)),'...Continued'), 1, 'last');
    events_details{last_valid_line} = [events_details{last_valid_line} events_details{continued_event_lines_IX(ii_continued_event)}];
end
events_IX(continued_event_lines_IX) = [];
events_TS(continued_event_lines_IX) = [];
events_TS_source(continued_event_lines_IX) = [];
events_type(continued_event_lines_IX) = [];
events_details(continued_event_lines_IX) = [];

repeating_IX = diff(events_IX)==0;
events_TS(repeating_IX) = [];
events_TS_source(repeating_IX) = [];
events_type(repeating_IX) = [];
events_details(repeating_IX) = [];

%% join event type with event details to a single string
EventStrings = {};
for ii_event = 1:length(events_type)
    event_str = [events_type{ii_event} ' - ' events_details{ii_event}];
    EventStrings{ii_event} = event_str;
end

%% extract recording parameteres from the event file

disp('2. reading+saving parameters from event file...');

is_recording = false;
for i=1:size(events_type,1)
    if strcmp(events_type{i},'Recording parameters')
        is_recording = true;
        paramst = strsplit(events_details{i},';');
        for i=1:length(paramst)
            if length(paramst{i})>1
                t = strsplit(paramst{i},' = ');
                
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
        break
    end
end


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
    ts_source_transceiver_events_IX = find(strcmp('Transceiver', events_TS_source)); % TODO: add Transceiver (fine) - validate string !!
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
    bv=str2num( curr_event_details(BV_values_str_interval) );
    if ~isempty(bv)
    BV_values(end+1) = bv;
    BV_timestamps(end+1) = events_TS(PC_gen_events_IX(ii_event));
    end
end
usec_2_min = 1/60*1e-6;
% BV_timestamps = BV_timestamps(4:end);
% BV_values = BV_values(4:end)
figure
plot( usec_2_min.*(BV_timestamps - BV_timestamps(1)), BV_values , '-o')
title('Battery discharge')
xlabel('time (minutes)')
ylabel('Voltage (V)')
grid on

hold on
if is_recording
    mode_changed_events_ts = events_TS(mode_changed_events_IX) - BV_timestamps(1);
    ylimits = get(gca, 'ylim');
    plot( repmat(usec_2_min .* mode_changed_events_ts,1,2)' , ylimits', '--m' )
end

legend({'battey voltage','record mode change event'})

fig_file = fullfile(Nlx_OutDir, 'Battery_discharge');
saveas(gcf, fig_file , 'jpg')
saveas(gcf, fig_file , 'fig')


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
    file_str = regexp(file_str,'\d*','Match');
    if length(file_str{1:end})<4
        file_name = [DATA_file_prefix '0' file_str{1:end} '.DT4'];
    else
        file_name = [DATA_file_prefix file_str{1:end} '.DT4'];
    end
    
    fid = fopen(fullfile(Nlg_InDir, file_name));
    filedata = fread(fid, 'uint16', 0, 'l');
    fclose(fid);
    filedata = double(filedata);
    data = reshape(filedata, num_channels, length(filedata)/num_channels);
    file_TS_usec = FileStarted_TS(ii_file_start_entry);
    
    %% remove DC
    if is_remove_DC
        for ii_channel = 1:size(data,1)
            %         data(ii_channel,:) = data(ii_channel,:) - mean(data(ii_channel,:));  % TAMIR: we prefer not to use this since data could be biased around real zero voltage (artifacts for example...)
            data(ii_channel,:) = data(ii_channel,:) - zero_DC_level_bit;         % There should be a theoretical constant representing the real zero voltage, but because the INTAN have a DC dependency on freq it is not perfect
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
        file_name = ['CSC_TT' num2str(TT_name(cnl+1)) '_' num2str(ch_per_TT_name(cnl+1)) '.ncs'];
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

%%
disp('Finished converting all files from Nlg format to Nlx format!')









%%