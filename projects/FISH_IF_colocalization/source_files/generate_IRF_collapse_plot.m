function generate_IRF_collapse_plot(FISH_IRF,IF_IRF,ind_flag,figure_name)
% Generates IRF plots
%   Function takes in FISH, IF IRF data inputs, as well as flag to plot
%   individual IRF lines. The figure is saved with figure_name
figure;
set(gcf, 'Visible', 'on');

subplot(2,1,1);
dist =[];
mean_FISH = [];
mean_IF =[];
for k =1:1:length(FISH_IRF)
    if ind_flag
        plot(smooth(FISH_IRF{k}.dist,15),smooth(FISH_IRF{k}.intensity,15),'LineWidth',0.25,'Color',[0.25 0 0]) ; hold on;
    end
    mean_FISH = [mean_FISH smooth(FISH_IRF{k}.intensity,15)]; 

end
mean_FISH = mean(mean_FISH,2);
plot(smooth(FISH_IRF{k}.dist,15), mean_FISH,'LineWidth',4,'Color',[1 0 0]);
subplot(2,1,2);
for k =1:1:length(FISH_IRF)
    if ind_flag
        plot(smooth(FISH_IRF{k}.dist,15),smooth(IF_IRF{k}.intensity,15),'LineWidth',0.25,'Color',[0 0.25 0]) ; hold on;
    end
    mean_IF = [mean_IF smooth(IF_IRF{k}.intensity,15)]; 
    
end
mean_IF = mean(mean_IF,2);
plot(smooth(FISH_IRF{k}.dist,15), mean_IF,'LineWidth',4,'Color',[0 1 0]);

set(gcf, 'Position', get(0, 'Screensize'));

C1 =corrcoef(mean_FISH,mean_IF);
disp(['Correlation coefficient for FISH and IF is' num2str(C1(1,2))]);    

% saveas(gcf,[figure_name '.fig']);
% saveas(gcf,[figure_name '.svg']);

end