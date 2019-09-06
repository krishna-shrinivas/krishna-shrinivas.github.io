function [Rand_protein,Rand_IF] = generate_average_images(Image_store,Fish_store,figure_name,IF_name,FISH_name)
% Generates average projection FISH/IF plots
%   Function takes in FISH, IF IRF image matrices, as well as data names.
%   The figure is saved under figure_name
%   Rand_protein and Rand_IF are just randomized permutations of the input
%   matrix that are one possible control

close all;
IF_proc = [];
FISH_proc = [];
for i=1:1:length(Fish_store)
    IF_proc(:,:,i) = mean(Image_store{i},3);
    FISH_proc(:,:,i) = mean(Fish_store{i},3);
end
average_FISH = mean(FISH_proc,3);
average_protein = mean(IF_proc,3);

for j=1:1:size(IF_proc,3)
   Temp_image = IF_proc(:,:,j);
   Rand_protein(:,:,j) = reshape(Temp_image(randperm(numel(Temp_image))),length(Temp_image),length(Temp_image));
   Temp_image = FISH_proc(:,:,j);
   Rand_IF(:,:,j) = reshape(Temp_image(randperm(numel(Temp_image))),length(Temp_image),length(Temp_image));

end
Rand_mean = mean(Rand_protein,3);
protein_color_map = [0 0.05 0;0 0.1 0; 0 0.15 0; 0 0.2 0;0 0.5 0;0 0.6 0;0 0.7 0;0 0.8 0;];
% FISH_color_map = [0.05 0 0;0.1 0 0; 0.15 0 0 ; 0.2 0 0 ;0.5 0 0 ;0.6 0 0 ;0.7 0 0 ;0.8 0 0 ;];

n_colors = 16;
FISH_color_map = zeros(n_colors,3);
FISH_color_map(:,1) = linspace(0,1,n_colors);


n_colors = 6;

protein_color_map = zeros(n_colors,3);
protein_color_map(:,2) = linspace(0,1,n_colors);


figure;
set(gcf, 'Visible', 'on');


ax1= subplot(2,2,1);

[Cf,hf]=contourf(average_FISH);
colormap(ax1,FISH_color_map); hold on;
hf.LineStyle='none';
% hf.LevelStep =30;
title([FISH_name ' FISH'],'FontSize',16);
pbaspect([1 1 1]);

ax2 = subplot(2,2,2);
surf(average_FISH); hold on;
view(ax2,[298.1 26]);
colormap(ax2,FISH_color_map);
title([FISH_name ' FISH intensity'],'FontSize',16);
colorbar
pbaspect([3 3 2]);

% colormap default
ax3=subplot(2,2,3);
[C,h]=contourf(average_protein); hold on;
h.LineStyle ='none';
% h.LevelStep =15;
colormap(ax3,protein_color_map); 
title([IF_name ' IF'],'FontSize',16);
pbaspect([1 1 1]);

ax4=subplot(2,2,4);
surf(average_protein); hold on;
view(ax4,[298.1 26]);
title([IF_name ' IF Intensity'],'FontSize',16);
colormap(ax4,protein_color_map); 
colorbar
pbaspect([3 3 2]);


% ax5=subplot(3,2,5);
% [Cr,hr]=contourf(Rand_mean); hold on;
% hr.LevelList = h.LevelList;
% hr.LineStyle ='none';
% colormap(ax5,protein_color_map); 
% title('Randomized IF','FontSize',16);
% 
% % hr.LevelStep =20;
% 
% 
% ax6 =subplot(3,2,6);
% surf(Rand_mean); hold on;
% colormap(protein_color_map);
% view(ax6,[298.1 26]);
% title('Randomized IF intensity','FontSize',16);
% colorbar

set(gcf, 'Position', get(0, 'Screensize'));
% 
% saveas(gcf,[figure_name '.fig']);
% saveas(gcf,[figure_name '.svg']);
end