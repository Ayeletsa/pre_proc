function nlg_TTL_ts_us=initial_extract_nlg_TTL(nlx_dir)

nlg_TTL_file_name = fullfile(nlx_dir, 'EVENTS__Digital in.nev');
FieldSelection = [1 0 0 0 0];
ExtractHeader = 0;
ExtractMode = 1;
ModeArray = [];
nlg_TTL_ts_us = Nlx2MatEV( nlg_TTL_file_name ,FieldSelection,ExtractHeader,ExtractMode,ModeArray);
end