function [L_Ratio,Isolation_dis] = cluster_quality(all_spikes_samples,all_spikes_ts,cluster_spikes_ts)

ind = [];
for n=1:length(cluster_spikes_ts)
    ind = [ind find(all_spikes_ts == cluster_spikes_ts(n))];
end
num_spikes_cluster = length(cluster_spikes_ts);
num_spikes_noise = length(all_spikes_ts) - num_spikes_cluster;


%% get feature space
all_spikes_samples = permute(all_spikes_samples,[3,1,2]); %Mx32x4
[nSpikes, nSamp, nch] = size(all_spikes_samples);

nrg = sqrt(sum(all_spikes_samples.^2,2))/nSamp; %energy per spikes per channel
Nrg_factor = repmat(nrg,1,32,1);
all_data_normalized = all_spikes_samples./Nrg_factor;

nrg = reshape(nrg,nSpikes,4); %energy per spikes per channel Mx4
for i = 1:4
    [coeffC,scores,latent] = pca(all_data_normalized(:,:,i));
    pc1(:,i) = scores(:,1);
end

fetures_space = [nrg pc1]; %Mx8
noise_feature_vect = fetures_space; noise_feature_vect(ind,:) = [];
cluster_feature_vect = fetures_space(ind,:);

%% Mahalanobis distance:
D_2_noise = mahal(noise_feature_vect,cluster_feature_vect);

%% L-Ratio:
CDF_noiseSpikes = cdf('Chisquare',D_2_noise,8); % 8 degrees of freedom
L_Ratio = sum(ones([num_spikes_noise,1]) - CDF_noiseSpikes)/num_spikes_cluster;

%% Isolation Distance
D_2_noise_sorted = sort(D_2_noise,'ascend');
if num_spikes_cluster<=length(D_2_noise_sorted)
    Isolation_dis = D_2_noise_sorted(num_spikes_cluster);
else
    Isolation_dis = D_2_noise_sorted(end);
end


