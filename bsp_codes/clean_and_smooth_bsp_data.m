function bsp_data = clean_and_smooth_bsp_data(bsp_data)

% TODO - measure balls distance
n_tags = length(bsp_data);
balls_pos = [0 125];

%% clear bsp samples (x,y,z) where -
% a. x is not between the balls
% b. y is far outside the hamama

half_hamama_width_m = 1.3;
y_std_criteria = 3;
x_dis_criteria = 10;
cleaning_ind = cell(1,n_tags);

for ii_tag = 1:n_tags
    
    x_cleaning_ind = or(bsp_data(ii_tag).normalized_pos(:,1) < (balls_pos(1) - x_dis_criteria), ...
    bsp_data(ii_tag).normalized_pos(:,1) > (balls_pos(2) + x_dis_criteria));
    
    y_std = std(bsp_data(ii_tag).normalized_pos(:,2));
    y_cleaning_criteria = half_hamama_width_m + y_std_criteria*y_std;
    y_cleaning_ind = (abs(bsp_data(ii_tag).normalized_pos(:,2)) > y_cleaning_criteria);
    
    cleaning_ind{ii_tag} = or(x_cleaning_ind,y_cleaning_ind);
    bsp_data(ii_tag).normalized_pos(cleaning_ind{ii_tag},:) = nan;
    
end


%% smooth the data
% 
% span = 3;
% u = [1 1 1]./3;
% for ii_tag = 1:n_tags
%      smoothed_x_pos =nan*zeros(size(bsp_data(ii_tag).normalized_pos(:,1)));
%     smoothed_y_pos =nan*zeros(size(bsp_data(ii_tag).normalized_pos(:,1)));
%     smoothed_z_pos = nan*zeros(size(bsp_data(ii_tag).normalized_pos(:,1)));
% 
%     x_pos = bsp_data(ii_tag).normalized_pos(:,1);
%     y_pos = bsp_data(ii_tag).normalized_pos(:,2);
%     z_pos = bsp_data(ii_tag).normalized_pos(:,3);
%      non_nan_ind=find(~isnan(x_pos));
%    % smoothed_x_pos = smooth(x_pos,span);
%    % smoothed_y_pos = smooth(y_pos,span);
%   %  smoothed_z_pos = smooth(z_pos,span);
%   
%     smoothed_x_pos(non_nan_ind) =conv(x_pos(non_nan_ind),u,'same');
%     smoothed_y_pos(non_nan_ind) = conv(y_pos(non_nan_ind),u,'same');
%     smoothed_z_pos(non_nan_ind) = conv(z_pos(non_nan_ind),u,'same');
% 
% 
%     
%     smoothed_pos = [smoothed_x_pos , smoothed_y_pos , smoothed_z_pos];
%     smoothed_pos(cleaning_ind{ii_tag},:) = nan;  % set again unwanted samples to NaNs after smoothing
%     bsp_data(ii_tag).normalized_pos = smoothed_pos;
%     
% end

end