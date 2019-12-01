function behavioral_modes=plot_behavior_and_correct_modes(behav_modes_plot,behavioral_modes,bsp_proc_data,tag_i,bat,day,behave_analysis_fig_dir_out,behav_params_file_name,dir_param_file_name,general_behavior_data_file_name,behave_day_struct_folder)
load(behav_params_file_name)
load(dir_param_file_name)
load(general_behavior_data_file_name)

us_factor=1e6;
number_of_plots=9;
%relevant times for plot:
FE_ind=bsp_proc_data(tag_i).flight_ind;
ts=bsp_proc_data(tag_i).ts;

total_time=ts(FE_ind(end))-ts(FE_ind(1));
time_bin=total_time/number_of_plots;
time_vec_for_plot=ts(FE_ind(1)):time_bin:ts(FE_ind(end));


%% figure param:
panel_size=[0.92 0.095];
x_position=0.02;
direc_color=[0 0.447 0.741;0.8500    0.3250    0.0980];

vertical_dist=0.1;
y_position(1)=0.85;
y_position(2)=y_position(1)-vertical_dist;
y_position(3)=y_position(2)-vertical_dist;
y_position(4)=y_position(3)-vertical_dist;
y_position(5)=y_position(4)-vertical_dist;
y_position(6)=y_position(5)-vertical_dist;
y_position(7)=y_position(6)-vertical_dist;
y_position(8)=y_position(7)-vertical_dist;
y_position(9)=y_position(8)-vertical_dist;

%% PLOT
figure('units','normalized','outerposition',[0 0 1 1])

% if tag_i==2
color_self=[0 1 0];
color_other=[1 0 0];

% title:
axes('position', [0.4, y_position(1), 0.1,  0.1]);
title(sprintf('Bat %d day %d - Total recording time %.1f min, Total flight time %.1f min',bat,day,(total_time/us_factor)/60,(length(FE_ind)/frame_per_second)/60), 'FontSize', 11, 'FontWeight', 'Bold');
box off;
axis off;
% else
%     color_self=[1 0 0];
%     color_other=[0 1 0];
% end
% All behavior
%--------------------------------------------------
behav_color=[0 0 1;0 1 0;1 0 0;0.5 0.5 0.5;1 1 0];

for time_i=1:length(time_vec_for_plot)-1
    ax(time_i)=axes('position',[x_position y_position(time_i) panel_size]);
    relevant_behav_ind=find(ts>=time_vec_for_plot(time_i) & ts<time_vec_for_plot(time_i+1));
    by_pass=behavioral_modes.bypass_ind(find(ts(behavioral_modes.bypass_ind)>=time_vec_for_plot(time_i) & ts(behavioral_modes.bypass_ind)<time_vec_for_plot(time_i+1)));
    plot(ts(relevant_behav_ind),pos_other_x(relevant_behav_ind),'.','color',color_other)
    hold on;
    plot(ts(relevant_behav_ind),pos_self_x(relevant_behav_ind),'.','color',color_self);
    hold on;
    if sum(isnan(pos_self_x(relevant_behav_ind)))==length(pos_self_x(relevant_behav_ind))
        set(gca,'xlim',[min(ts(relevant_behav_ind)) max(ts(relevant_behav_ind))])
    else
        set(gca,'xlim',[min(ts(relevant_behav_ind)) max(ts(relevant_behav_ind))],'ylim' ,[min(pos_self_x(relevant_behav_ind)) max(pos_self_x(relevant_behav_ind))])
    end
    set(gca,'XTick',[min(ts) max(ts)],'XTickLabel',round(([min(ts) max(ts)])/1e3/60))
    if ~isempty(by_pass)
        plot(ts(by_pass),pos_self_x(by_pass),'*k');
    end
    
    for behav_mod_i=1:length(behav_modes_plot)
        for ii=1:length(behav_modes_plot(behav_mod_i).start)
            start_point=[];
            end_point=[];
            if ts(behav_modes_plot(behav_mod_i).start(ii))>=time_vec_for_plot(time_i) & ts(behav_modes_plot(behav_mod_i).end(ii))<=time_vec_for_plot(time_i+1)
                start_point=ts(behav_modes_plot(behav_mod_i).start(ii));
                end_point=ts(behav_modes_plot(behav_mod_i).end(ii));
            elseif ts(behav_modes_plot(behav_mod_i).start(ii))>=time_vec_for_plot(time_i) & ts(behav_modes_plot(behav_mod_i).end(ii))>=time_vec_for_plot(time_i+1)
                start_point=ts(behav_modes_plot(behav_mod_i).start(ii));
                end_point=time_vec_for_plot(time_i+1);
            elseif ts(behav_modes_plot(behav_mod_i).start(ii))<=time_vec_for_plot(time_i) & ts(behav_modes_plot(behav_mod_i).end(ii))<=time_vec_for_plot(time_i+1)
                start_point=time_vec_for_plot(time_i);
                end_point=ts(behav_modes_plot(behav_mod_i).end(ii));
                
                
            end
            if ~isempty(start_point)
                rec_x=[start_point  end_point  end_point start_point];
                rec_y=[min(pos_self_x)  min(pos_self_x) max(pos_self_x) max(pos_self_x)];
                p=patch(rec_x,rec_y,behav_color(behav_mod_i,:),'EdgeColor','none');
                set(p,'FaceAlpha',0.3)
            end
        end
    end
end

%% manual corrections
if correct_manually
    
    % 1. Add CO
    q1=str2num(cell2mat(inputdlg('Do you want to add CO? (1/0)')));
    if q1
        how_many_CO=str2num(cell2mat(inputdlg('How many CO do you want to add?')));
        [co_time_to_add,~]=ginput(how_many_CO);
        [~,ind]=min(abs(ts-co_time_to_add'));
        [val,indx]=min(abs(distance_change_sign-ind));
        co_point_to_add=distance_change_sign(indx(val<manual_min_dis_from_CO)); %add CO only if it found closest point of changed sign
        behavioral_modes.CO_point=sort([behavioral_modes.CO_point co_point_to_add']);
        ind_to_add=[];
        co_event={};
        for co_i=1:length(co_point_to_add)
            co=co_point_to_add(co_i);
            %find ax:
            poss_ax=find([time_vec_for_plot]<ts(co));
            axs=poss_ax(end);
            relevant_ax=ax(axs);
            axes(relevant_ax)
            plot(ts(co),pos_self_x(co),'ko')
            
            %find window:
            % find the start:
            long_dist_ind=find(abs(distnace_other_from_self)>=CO_window(1));
            before_CO=long_dist_ind(find([long_dist_ind-co]<0));
            start_ind=before_CO(end);%first ind before CO that the distance between the bats is lrager than ..
            if (co-start_ind)>max_wind_CO
                start_ind=co-max_wind_CO;
            end
            %find the end of the window:
            long_dist_ind=find(abs(distnace_other_from_self)>=CO_window(2));
            after_CO=long_dist_ind(find([long_dist_ind-co]>0));
            if ~isempty(after_CO)
                end_ind=after_CO(1); %first ind after CO that the distance between the bats is lrager than ..
                
                behavioral_modes.CO_ind=sort([behavioral_modes.CO_ind start_ind:end_ind]);
                behavioral_modes.CO_event{end+1}=start_ind:end_ind;
            end
            
        end
        manually_added.co_point_to_add=co_point_to_add;
    else
        manually_added.co_point_to_add=[];
    end
    
    %2. Add solo
    q2=str2num(cell2mat(inputdlg('Do you want to add solo? (1/0)')));
    if q2
        how_many_solo=str2num(cell2mat(inputdlg('How many Solo do you want to add?')));
        manually_added.solo_ind_to_add=[];
        for solo_i=1:how_many_solo
            [solo_time_to_add,~]=ginput(2);
            [~,ind]=min(abs(ts-solo_time_to_add'));
            relevant_ind_to_add=ind(1):ind(2);
            relevant_ind_to_add=setdiff(relevant_ind_to_add,behavioral_modes.solo_ind);
            behavioral_modes.solo_ind=sort([behavioral_modes.solo_ind relevant_ind_to_add]);
            manually_added.solo_ind_to_add=[manually_added.solo_ind_to_add,relevant_ind_to_add];
            %find ax:
            poss_ax=find([time_vec_for_plot]<ts(relevant_ind_to_add(1)));
            axs=poss_ax(end);
            relevant_ax=ax(axs);
            axes(relevant_ax)
            plot([ts(relevant_ind_to_add(1)) ts(relevant_ind_to_add(end))],[ pos_self_x(relevant_ind_to_add(1)) pos_self_x(relevant_ind_to_add(end))],'k+')
        end
    else
        manually_added.solo_ind_to_add=[];
    end
    
    % 3 remove solo
    q3=str2num(cell2mat(inputdlg('Do you want to remove solo? (1/0)')));
    if q3
        how_many_solo=str2num(cell2mat(inputdlg('How many Solo do you want to remove?')));
        manually_added.solo_ind_to_remove=[];
        for solo_i=1:how_many_solo
            [solo_time_to_add,~]=ginput(2);
            [~,ind]=min(abs(ts-solo_time_to_add'));
            relevant_ind_to_remove=[min(ind) max(ind)];
            solo_ind_cleaned=setdiff(behavioral_modes.solo_ind,relevant_ind_to_remove);
            behavioral_modes.solo_ind=solo_ind_cleaned;
            manually_added.solo_ind_removed=relevant_ind_to_remove;
            %find ax:
            poss_ax=find([time_vec_for_plot]<ts(relevant_ind_to_remove(1)));
            axs=poss_ax(end);
            relevant_ax=ax(axs);
            axes(relevant_ax)
            plot([ts(relevant_ind_to_remove(1)) ts(relevant_ind_to_remove(end))],[ pos_self_x(relevant_ind_to_remove(1)) pos_self_x(relevant_ind_to_remove(end))],'k+')
        end
    else
        manually_added.solo_ind_removed=[];
    end
    
    file_name=fullfile(behave_day_struct_folder,['manually_added_behavioral_modes_bat_',num2str(bat),'_day_',num2str(day),'.mat']);
    save(file_name,'manually_added')
else
    file_name=fullfile(behave_day_struct_folder,['manually_added_behavioral_modes_bat_',num2str(bat),'_day_',num2str(day),'.mat']);
    load(file_name)
%     %add:
%     relev_fields=fields(manually_added);
%     
%     if strcmp(relev_fields,'solo_ind_to_add')
%         behavioral_modes.solo_ind=sort([behavioral_modes.solo_ind manually_added.solo_ind_to_add]);
%     end
%     if strcmp(relev_fields,'solo_ind_removed')
%         behavioral_modes.solo_ind=setdiff(behavioral_modes.solo_ind,manually_added.solo_ind_removed);
%         
%     end
%     if strcmp(relev_fields,'co_ind_to_add')
%         behavioral_modes.CO_point=sort([behavioral_modes.CO_point manually_added.co_ind_to_add]);
%     end
end
%% create time scale:
axes('position',[x_position 0.03 panel_size(1) 0.02]);
time_bin_in_min=(time_bin/us_factor)/60;
plot([0 time_bin/time_bin_in_min], [0 0],'k','LineWidth',5)
xlim([0 time_bin])
text(time_bin/(2*time_bin_in_min),-1.5, '1 min','HorizontalAlignment','center')
box off
axis off
%create legend
axes('position',[0.7 0.95 0.2 0.05])
x_pos=0:1/length(behav_modes_plot):1;
rec_size=0.8*(1/length(behav_modes_plot));
for behav_mod_i=1:length(behav_modes_plot)
    rec_x=[x_pos(behav_mod_i) x_pos(behav_mod_i)+rec_size x_pos(behav_mod_i)+rec_size x_pos(behav_mod_i)];
    rec_y=[0  0 1 1];
    p=patch(rec_x,rec_y,behav_color(behav_mod_i,:),'EdgeColor','none');
    set(p,'FaceAlpha',0.3)
    if behav_mod_i==3
        text(rec_x(1)+0.01,0.7,behav_modes_plot(behav_mod_i).name(1:5))
        text(rec_x(1)+0.01,0.3,behav_modes_plot(behav_mod_i).name(7:13))
    else
        text(rec_x(1)+0.01,0.5,behav_modes_plot(behav_mod_i).name)
    end
    box off
    axis off
end

% save
fig_name=fullfile(behave_analysis_fig_dir_out,['behavioral_mode_bat_',num2str(bat),'_day_',num2str(day),'.tif']);
saveas(gcf,fig_name)
fig_name=fullfile(behave_analysis_fig_dir_out,['behavioral_mode_bat_',num2str(bat),'_day_',num2str(day),'.fig']);
saveas(gcf,fig_name)
clf
