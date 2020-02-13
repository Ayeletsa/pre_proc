function bsp_TTL_ts_ns=initial_extract_bsp_TTL(bsp_dir)

%1. extract bsp data convert files
bsp_extract_data (bsp_dir);

%2. readr TTLs:
load( fullfile(bsp_dir, 'bsp_TTL.mat') );
bsp_TTL_ts_ns = bsp_TTL_ts_ns';
end



