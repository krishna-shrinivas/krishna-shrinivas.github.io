function production_quality_IRF(FISH_IRF,IF_IRF,Random_IRF,ind_flag,figure_name,FISH_name,IF_name,FISH_color,IF_color)
%   Produce "production" quality IRF plots from input data
%   This function was broadly tweaked to generate the plots that were
%   eventually used for publication. The inputs are the IRF functions for
%   the FISH and IF, as well as random IF data. Various labels are input
%   parameters, and the color choice.
%   Figures are saved under figure_name
input_params.xpixel = 0.0572;

input_params.zpixel = 0.2;

% for i =1:1:length(FISH_data)
%     [FISH_IRF{i}] = convert_3D_image_IRF(FISH_data{i},input_params);
%     [IF_IRF{i}] = convert_3D_image_IRF(IF_data{i},input_params);
% end
figure1=figure;
set(gcf, 'Visible', 'on');
green_color = [0 0.800000011920929 0.400000005960464];
magenta_color = [1 0 1];
red_color = [1 0 0];
set(figure1,'defaultAxesColorOrder',[FISH_color; IF_color]);

subplot11 = subplot(1,1,1,'Parent',figure1);
set(subplot11,'defaultAxesColorOrder',[FISH_color; IF_color]);

hold(subplot11,'on');
dist =[];
mean_FISH = [];
mean_IF =[];
yyaxis left

for k =1:1:length(FISH_IRF)
    if ind_flag
        plot(smooth(FISH_IRF{k}.dist,15),smooth(FISH_IRF{k}.intensity,15),'LineWidth',0.25,'Color',FISH_color*0.5) ; hold on;
    end
    mean_FISH = [mean_FISH smooth(FISH_IRF{k}.intensity,15)]; 

end
mean_FISH = mean(mean_FISH,2);
plot(smooth(FISH_IRF{k}.dist,15), mean_FISH,'LineWidth',4,'Color',FISH_color);
ylabel([FISH_name ' FISH']);

yyaxis right
for k =1:1:length(FISH_IRF)
    if ind_flag
        plot(smooth(FISH_IRF{k}.dist,15),smooth(IF_IRF{k}.intensity,15),'LineWidth',0.25,'Color',IF_color*0.5) ; hold on;
    end
    mean_IF = [mean_IF smooth(IF_IRF{k}.intensity,15)]; 
    
end
mean_IF = mean(mean_IF,2);
p1=plot(smooth(FISH_IRF{k}.dist,15), mean_IF,'LineWidth',4,'Color',IF_color);
ylabel([IF_name ' IF']);
set(subplot11,'FontSize',24);

if ~isempty(Random_IRF)
    mean_R_IF = [];
    yyaxis right
    for k =1:1:length(Random_IRF)
        if ind_flag
            plot(smooth(Random_IRF{k}.dist,15),smooth(Random_IRF{k}.intensity,15),'LineWidth',0.25,'Color',IF_color*0.5,'LineStyle','--') ; hold on;
        end
        mean_R_IF = [mean_R_IF smooth(Random_IRF{k}.intensity,15)]; 

    end
    mean_R_IF = mean(mean_R_IF,2);
    p2=plot(smooth(Random_IRF{k}.dist,15), mean_R_IF,'LineWidth',4,'Color',IF_color*0.4,'LineStyle','--');
    C1 =corrcoef(mean_FISH,mean_IF);

    legend([p1,p2],{['\rho = ' num2str(C1(1,2),2) ' N_{foci} = ' num2str(length(FISH_IRF))],[IF_name 'Random IF']});
    set(subplot11,'FontSize',24);
else
    C1 =corrcoef(mean_FISH,mean_IF);

  legend([p1],{['\rho = ' num2str(C1(1,2),2) ' N_{foci} = ' num2str(length(FISH_IRF))]});

end
set(gcf, 'Position', get(0, 'Screensize'));
xlabel('\mum');

saveas(gcf,[figure_name '.fig']);
saveas(gcf,[figure_name '.svg']);

% 
% 
% 
% 
% 













end

function old_IRF(N_Image_store,N_Fish_store,N_random_IF,idf_name,IF_name,FISH_name,bg_subtract,max_multiplier)

% close all;


pixel_to_micron = 1/0.0572;


%% Generate data and plot for each trajectory

center_point = [ceil(size(N_Image_store,2)/2),ceil(size(N_Image_store,2)/2)];
dist_FISH = [];
count =1;
temp_store_FISH = [];
temp_store_IF = [];


temp_store_normalized_FISH = [];
temp_store_normalized_IF = [];



figure1=figure;
% set(gcf, 'Visible', 'off');

%% Generate for the average plots
count =1;
if ~bg_subtract
    average_protein = mean(N_Image_store,3);
    average_FISH = mean(N_Fish_store,3);
    average_random_protein = mean(N_random_IF,3);
else
    average_FISH = mean(N_Fish_store,3);
    average_protein = mean(N_Image_store,3) - mean(mean(mean(N_random_IF,3)));
    average_random_protein = mean(N_random_IF,3)- mean(mean(mean(N_random_IF,3)));

end
dist_max = 0.5*max_multiplier*size(average_FISH,2)/pixel_to_micron;

for i =1:1:size(average_FISH,1)
    for j=1:1:size(average_FISH,2)
        dist_FISH(count)= sqrt((i-center_point(1))^2+(j-center_point(2))^2)/pixel_to_micron;
        if dist_FISH(count)<dist_max
            Intensity_FISH(count) = average_FISH(i,j);
            Intensity_protein(count) = average_protein(i,j);
            Intensity_random_protein(count) = average_random_protein(i,j);
            count = count+1;
        else
            dist_FISH(count) = [];        
        end
       
    end
end
[sort_dist_FISH,index_sort]=sort(dist_FISH);
[sort_Intensity_FISH] = Intensity_FISH(index_sort);
[sort_Intensity_protein] = Intensity_protein(index_sort);
[sort_Intensity_random_protein] = Intensity_random_protein(index_sort);


Normalized_by_mean_FISH = sort_Intensity_FISH/ mean(sort_Intensity_FISH);
Normalized_by_mean_protein = sort_Intensity_protein/ mean(sort_Intensity_protein);
Normalized_by_mean_random_protein = sort_Intensity_random_protein/ mean(sort_Intensity_protein);

d= 1;
while(d <= length(sort_dist_FISH))
    repeats = find(sort_dist_FISH== sort_dist_FISH(d));
    if length(repeats)~=1
        sort_dist_FISH(repeats) = [];
        sort_Intensity_FISH(d) = mean(sort_Intensity_FISH([d, repeats]));
        Normalized_by_mean_FISH(d) = mean(Normalized_by_mean_FISH([d, repeats]));

        
        sort_Intensity_protein(d) = mean(sort_Intensity_protein([d, repeats]));
        Normalized_by_mean_protein(d) = mean(Normalized_by_mean_protein([d, repeats]));

        sort_Intensity_random_protein(d) = mean(sort_Intensity_random_protein([d, repeats]));
        Normalized_by_mean_random_protein(d) = mean(Normalized_by_mean_random_protein([d, repeats]));

        
        sort_Intensity_FISH(repeats) = [];
        Normalized_by_mean_FISH(repeats) = [];
       
        
        sort_Intensity_protein(repeats) = [];
        Normalized_by_mean_protein(repeats) = [];
       
        sort_Intensity_random_protein(repeats) = [];
        Normalized_by_mean_random_protein(repeats) = [];


    end
    d=d+1;
end

C1 =corrcoef(sort_Intensity_FISH,sort_Intensity_protein);
C2 =corrcoef(sort_Intensity_FISH,sort_Intensity_random_protein);
disp(['Correlation coefficient for FISH and IF is' num2str(C1(1,2))]);    
disp(['Correlation coefficient for FISH and random IF is' num2str(C2(1,2))]);    

green_color = [0 0.800000011920929 0.400000005960464];
magenta_color = [1 0 1];
set(figure1,'defaultAxesColorOrder',[magenta_color; green_color]);

subplot11 = subplot(2,1,1,'Parent',figure1);
set(subplot11,'defaultAxesColorOrder',[magenta_color; green_color]);

hold(subplot11,'on');

% Create multiple lines using matrix input to plot
yyaxis left
plot(smooth(sort_dist_FISH,5),smooth(sort_Intensity_FISH,5),'LineWidth',4,'Color',magenta_color); hold on;
ylim([min(smooth(sort_Intensity_FISH,5)) max(smooth(sort_Intensity_FISH,5))]);
ylabel([FISH_name ' FISH']);

yyaxis right
plot(smooth(sort_dist_FISH,5),smooth(sort_Intensity_protein,5),'LineWidth',4,'Color',green_color); hold on;
plot(smooth(sort_dist_FISH,5),smooth(sort_Intensity_random_protein,5),'LineWidth',4,'LineStyle',':','Color',green_color); hold on;
ylim([min(min(smooth(sort_Intensity_random_protein,5)),min(smooth(sort_Intensity_protein,5))) max(max(smooth(sort_Intensity_protein,5)),max(smooth(sort_Intensity_random_protein,5)))]);
ylabel([IF_name ' IF']);

% Create xlabel
xlabel('\mum');
xlim([0 dist_max]);

% Create ylabel
%     ylim([0 1])
box(subplot11,'on');
% Set the remaining axes properties
set(subplot11,'FontSize',24);

subplot12 = subplot(2,1,2,'Parent',figure1);
set(subplot12,'defaultAxesColorOrder',[magenta_color; green_color]);

hold(subplot12,'on');
pixel_to_micron = 1/0.0572;

% Create multiple lines using matrix input to plot
yyaxis left

plot(smooth(sort_dist_FISH,5),smooth(Normalized_by_mean_FISH,5),'LineWidth',4,'Color',magenta_color); hold on;
ylim([min(smooth(Normalized_by_mean_FISH,5)) max(smooth(Normalized_by_mean_FISH,5))]);
ylabel([FISH_name ' FISH (Normalized)']);

yyaxis right
plot(smooth(sort_dist_FISH,5),smooth(Normalized_by_mean_protein,5),'LineWidth',4,'Color',green_color);
plot(smooth(sort_dist_FISH,5),smooth(Normalized_by_mean_random_protein,5),'LineWidth',4,'LineStyle',':','Color',green_color); hold on;
ylim([min(min(smooth(Normalized_by_mean_random_protein,5)),min(smooth(Normalized_by_mean_protein,5))) max(max(smooth(Normalized_by_mean_protein,5)),max(smooth(Normalized_by_mean_random_protein,5)))]);
ylabel([IF_name ' IF (Normalized)']);

% Create xlabel
xlabel('\mum');
xlim([0 dist_max]);

% Create ylabel
%     ylim([0 1])
box(subplot12,'on');
% Set the remaining axes properties
set(subplot12,'FontSize',24);




set(gcf, 'Position', get(0, 'Screensize'));

% saveas(gcf,[idf_name],'fig');
saveas(gcf,[idf_name '.svg'],'svg');

end