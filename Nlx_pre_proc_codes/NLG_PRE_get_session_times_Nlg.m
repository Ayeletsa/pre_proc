function p = NLG_PRE_get_session_times_Nlg(p)


%% In Nlg format - for each session - find the corresponding timestamps
% of its start and end, and update the structure. 
% Both Nlx and Nlg timestamps are in microseconds(!)

%----------------------------------------------------------------------------------------
for nses = 1:length(p.S)
    s = p.S(nses);
    
        start_time = p.Nlg_EventTimestamps(s.events(1))+s.time_offsets_in_seconds(1)*1e6;
        end_time = p.Nlg_EventTimestamps(s.events(2))+s.time_offsets_in_seconds(2)*1e6;
        p.S(nses).start_time = start_time;
        p.S(nses).end_time = end_time;
        
        if s.time_offsets_in_seconds(1) == 0
            start_event = p.Nlg_EventStrings{s.events(1)};
        elseif s.time_offsets_in_seconds(1) > 0
            start_event = sprintf('%s+%d',p.Nlg_EventStrings{s.events(1)},...
                s.time_offsets_in_seconds(1));
        else
            start_event = sprintf('%s-%d',p.Nlg_EventStrings{s.events(1)},...
                -s.time_offsets_in_seconds(1));
        end
        if s.time_offsets_in_seconds(2) == 0
            end_event = p.Nlg_EventStrings{s.events(2)};
        elseif s.time_offsets_in_seconds(2) > 0
            end_event = sprintf('%s+%d',p.Nlg_EventStrings{s.events(2)},...
                s.time_offsets_in_seconds(2));
        else
            end_event = sprintf('%s-%d',p.Nlg_EventStrings{s.events(2)},...
                -s.time_offsets_in_seconds(2));
        end
        
        
        p.S(nses).start_event = start_event;
        p.S(nses).end_event = end_event;
            
        % if we want to  use manually entered timestamp (in microseconds) (good for cases when we are missing the exact event)

        if s.exact_time_in_microseconds(1) ~=0
            p.S(nses).start_time = s.exact_time_in_microseconds(1);
            p.S(nses).start_event=sprintf('User defined %d (microseconds)',p.S(nses).start_time);
        end
        if s.exact_time_in_microseconds(2) ~=0
            p.S(nses).end_time =  s.exact_time_in_microseconds(2);
            p.S(nses).end_event=sprintf('User defined %d (microseconds)',p.S(nses).end_time);
        end
end



