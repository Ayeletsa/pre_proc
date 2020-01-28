s=fieldnames(p.bsp_data);

for ii=1:length(s)
    f=s{ii};
    
    for tag_i=1:2
        if ismember(ii,[1,10,11,18])
            new(tag_i).(f)=p.bsp_data(tag_i).(f);
            continue
        else
            x1=p.bsp_data(tag_i).(f);
            x2=p_b.bsp_data(tag_i).(f);
            if size(x2,1)>size(x2,2)
                new(tag_i).(f)=[x1;x2];
            else
                new(tag_i).(f)=[x1,x2];
                
            end
        end
        
    end
end