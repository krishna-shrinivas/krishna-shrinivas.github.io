function generate_foci_image(FISH,IF,figure_name,IF_name,FISH_name)
%   Generates contour plot of 3-D matrix with labels provided
close all;
protein_color_map = [0 0.05 0;0 0.1 0; 0 0.15 0; 0 0.2 0;0 0.5 0;0 0.6 0;0 0.7 0;0 0.8 0;];

n_colors = 16;
FISH_color_map = zeros(n_colors,3);
FISH_color_map(:,1) = linspace(0,1,n_colors);


n_colors = 6;

protein_color_map = zeros(n_colors,3);
protein_color_map(:,2) = linspace(0,1,n_colors);

figure;
set(gcf, 'Visible', 'off');

ax1= subplot(3,2,1);

[Cf,hf]=contourf(FISH);
colormap(ax1,FISH_color_map); hold on;
hf.LineStyle='none';
hf.LevelStep =30;
title([FISH_name ' FISH'],'FontSize',16);

ax2 = subplot(2,2,2);
surf(FISH); hold on;
view(ax2,[298.1 26]);
colormap(ax2,FISH_color_map);
title([FISH_name ' FISH intensity'],'FontSize',16);
colorbar

% colormap default
ax3=subplot(2,2,3);
[C,h]=contourf(IF); hold on;
h.LineStyle ='none';
h.LevelStep =15;
colormap(ax3,protein_color_map); 
title([IF_name ' IF'],'FontSize',16);

ax4=subplot(2,2,4);
surf(IF); hold on;
view(ax4,[298.1 26]);
title([IF_name ' IF Intensity'],'FontSize',16);
colormap(ax4,protein_color_map); 
colorbar


set(gcf, 'Position', get(0, 'Screensize'));

saveas(gcf,[figure_name '.fig']);
saveas(gcf,[figure_name '.svg']);
close all;

end