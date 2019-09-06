function  visualisation_data(filename)

load(filename)

%Averages cluster size over multiple trajectories for same parameters
for i=1:length(time)
    for p =1:length(time{1})
        %Tracks the mean cluster size for each trajectory
        %(i,p) represents (parameter choice, trajectory replicate)
        mean_cluster_size(i,p) = mean(cluster_max_size{i}{p});

    end
end


figure;
plot(Nc, mean(mean_cluster_size,2)./params.Nc,'LineWidth',4); hold on;
plot(Nc, (mean(mean_cluster_size,2)+2*std(mean_cluster_size')')./params.Nc,'--','LineWidth',2); hold on;
plot(Nc, (mean(mean_cluster_size,2)-2*std(mean_cluster_size')')./params.Nc,'--','LineWidth',2);hold on;

end