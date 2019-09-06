function production_average_images(Image_store,Fish_store,Random_store,figure_name,IF_name,FISH_name,FISH_color,IF_color,stack_range)
%   Produce "production" quality FISH-IF plots from input data
%   This function was broadly tweaked to generate the plots that were
%   eventually used for publication. The inputs are the image matrices for
%   the FISH and IF, as well as random IF data. Various labels are input
%   parameters, and the color choice.
%   Figures are saved under figure_name
pixeltoumxy = 0.0572;

if isempty(Random_store)
    close all;
    IF_proc = [];
    FISH_proc = [];
    for i=1:1:length(Fish_store)
        IF_proc(:,:,i) = mean(Image_store{i}(:,:,stack_range(1):stack_range(2)),3);
        FISH_proc(:,:,i) = mean(Fish_store{i}(:,:,stack_range(1):stack_range(2)),3);
    end
    average_FISH = mean(FISH_proc,3);
    average_protein = mean(IF_proc,3);

    n_colors = 16;
    FISH_color_map = zeros(n_colors,3);
    protein_color_map = zeros(n_colors,3);

    for i=1:1:3
        FISH_color_map(:,i) = linspace(0,FISH_color(i),n_colors);
        protein_color_map(:,i) = linspace(0,IF_color(i),n_colors);

    end

    % Color scheme 1
    hex_list = {'#000000','#18091f','#250d38','#330c52','#41066f','#480381','#3f0c7d','#351179','#2a1575','#1e1871','#37386d','#4c6765','#529859','#48ca44','#00ff00'};
    cmap_trial =hex2rgb(hex_list);

    figure;
    set(gcf, 'Visible', 'on');
    X = [0:1:size(average_protein,2)-1]*pixeltoumxy -0.5*size(average_protein,2)*pixeltoumxy;
    Y = [0:1:size(average_protein,2)-1]*pixeltoumxy -0.5*size(average_protein,2)*pixeltoumxy;

    ax1= subplot(2,2,1);

    [Cf,hf]=contourf(X,Y,average_FISH);
    colormap(ax1,FISH_color_map); hold on;
    xlabel('\mum');
    ylabel('\mum');
    pbaspect([2 2 1]);
    hf.LineStyle='none';
    % hf.LevelStep =30;
    title([FISH_name ' FISH'],'FontSize',16);

    minz = min(min(average_FISH));
    maxz = max(max(average_FISH));
    ax2 = subplot(2,2,2);
    surf(X,Y,average_FISH); hold on;
    xlabel('\mum');
    ylabel('\mum');

    view(ax2,[298.1 26]);
    zlim([minz maxz]);
    xlim([min(X) max(X)]);
    ylim([min(Y) max(Y)]);

    colormap(ax2,FISH_color_map);
    title([FISH_name ' FISH intensity'],'FontSize',16);
    colorbar

    % colormap default
    ax3=subplot(2,2,3);
    [C,h]=contourf(X,Y,average_protein); hold on;
    xlabel('\mum');
    ylabel('\mum');


    h.LineStyle ='none';
    % h.LevelStep =15;
    colormap(ax3,cmap_trial); 
    title([IF_name ' IF'],'FontSize',16);
    pbaspect([1 1 1]);

    minz = min(min(average_protein));
    maxz = max(max(average_protein));
    ax4=subplot(2,2,4);
    surf(X,Y,average_protein); hold on;
    view(ax4,[298.1 26]);
    xlabel('\mum');
    ylabel('\mum');

    zlim([minz maxz]);
    xlim([min(X) max(X)]);
    ylim([min(Y) max(Y)]);
    title([IF_name ' IF Intensity'],'FontSize',16);
    colormap(ax4,cmap_trial); 
    colorbar
else
    close all;
    IF_proc = [];
    FISH_proc = [];
    random_proc = []
    for i=1:1:length(Fish_store)
        IF_proc(:,:,i) = mean(Image_store{i}(:,:,stack_range(1):stack_range(2)),3);
        FISH_proc(:,:,i) = mean(Fish_store{i}(:,:,stack_range(1):stack_range(2)),3);
    end
    for i=1:1:length(Random_store)
        random_proc(:,:,i) = mean(Random_store{i}(:,:,stack_range(1):stack_range(2)),3);
    
    end
    average_FISH = mean(FISH_proc,3);
    average_protein = mean(IF_proc,3);
    average_random_protein = mean(random_proc,3);
    
    n_colors = 16;
    FISH_color_map = zeros(n_colors,3);
    protein_color_map = zeros(n_colors,3);

    for i=1:1:3
        FISH_color_map(:,i) = linspace(0,FISH_color(i),n_colors);
        protein_color_map(:,i) = linspace(0,IF_color(i),n_colors);

    end


    hex_list = {'#000000','#18091f','#250d38','#330c52','#41066f','#480381','#3f0c7d','#351179','#2a1575','#1e1871','#37386d','#4c6765','#529859','#48ca44','#00ff00'};
    cmap_trial =hex2rgb(hex_list);

    figure;
    set(gcf, 'Visible', 'on');
    X = [0:1:size(average_protein,2)-1]*pixeltoumxy -0.5*size(average_protein,2)*pixeltoumxy;
    Y = [0:1:size(average_protein,2)-1]*pixeltoumxy -0.5*size(average_protein,2)*pixeltoumxy;

    ax1= subplot(3,2,1);

    [Cf,hf]=contourf(X,Y,average_FISH);
    colormap(ax1,FISH_color_map); hold on;
    xlabel('\mum');
    ylabel('\mum');
    pbaspect([2 2 1]);
    hf.LineStyle='none';
    colorbar;
    % hf.LevelStep =30;
    title([FISH_name ' FISH'],'FontSize',16);

    minz = min(min(average_FISH));
    maxz = max(max(average_FISH));

    ax2 = subplot(3,2,2);
    surf(X,Y,average_FISH); hold on;
    xlabel('\mum');
    ylabel('\mum');

    view(ax2,[298.1 26]);
    zlim([minz maxz]);
    xlim([min(X) max(X)]);
    ylim([min(Y) max(Y)]);

    colormap(ax2,FISH_color_map);
    title([FISH_name ' FISH intensity'],'FontSize',16);
    colorbar

    
    minz = min(min(min(average_protein)),min(min(average_random_protein)));
    maxz = max(max(max(average_protein)),max(max(average_random_protein)));

    
    % colormap default
    ax3=subplot(3,2,3);
    [C,h]=contourf(X,Y,average_protein); hold on;
    xlabel('\mum');
    ylabel('\mum');

    h.LineStyle ='none';
    % h.LevelStep =15;
    ax3.CLim = [minz maxz];
    colormap(ax3,cmap_trial); 
    title([IF_name ' IF'],'FontSize',16);
    colorbar;

    pbaspect([1 1 1]);
   
    ax4=subplot(3,2,4);
    surf(X,Y,average_protein); hold on;
    view(ax4,[298.1 26]);
    xlabel('\mum');
    ylabel('\mum');
        ax4.CLim = [minz maxz];

    zlim([minz maxz]);
    xlim([min(X) max(X)]);
    ylim([min(Y) max(Y)]);
    title([IF_name ' IF Intensity'],'FontSize',16);
    colormap(ax4,cmap_trial); 
    colorbar
    
    ax5=subplot(3,2,5);
    [C,h]=contourf(X,Y,average_random_protein); hold on;
    xlabel('\mum');
    ylabel('\mum');
    ax5.CLim = [minz maxz];

    h.LineStyle ='none';
    % h.LevelStep =15;
    
    colormap(ax5,cmap_trial); 
    title([IF_name ' Random IF'],'FontSize',16);
    pbaspect([1 1 1]);
    colorbar;


    ax6=subplot(3,2,6);
    surf(X,Y,average_random_protein); hold on;
    view(ax6,[298.1 26]);
    xlabel('\mum');
    ylabel('\mum');

    zlim([minz maxz]);
    xlim([min(X) max(X)]);
    ylim([min(Y) max(Y)]);
        ax6.CLim = [minz maxz];

    title([IF_name ' Random IF Intensity'],'FontSize',16);
    colormap(ax6,cmap_trial); 
    colorbar

    
    
    
end

set(gcf, 'Position', get(0, 'Screensize'));

saveas(gcf,[figure_name '.fig']);
saveas(gcf,[figure_name '.svg']);
end

function [ rgb ] = hex2rgb(hex,range)
% hex2rgb converts hex color values to rgb arrays on the range 0 to 1. 
% 
% 
% * * * * * * * * * * * * * * * * * * * * 
% SYNTAX:
% rgb = hex2rgb(hex) returns rgb color values in an n x 3 array. Values are
%                    scaled from 0 to 1 by default. 
%                    
% rgb = hex2rgb(hex,256) returns RGB values scaled from 0 to 255. 
% 
% 
% * * * * * * * * * * * * * * * * * * * * 
% EXAMPLES: 
% 
% myrgbvalue = hex2rgb('#334D66')
%    = 0.2000    0.3020    0.4000
% 
% 
% myrgbvalue = hex2rgb('334D66')  % <-the # sign is optional 
%    = 0.2000    0.3020    0.4000
% 
%
% myRGBvalue = hex2rgb('#334D66',256)
%    = 51    77   102
% 
% 
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myrgbvalues = hex2rgb(myhexvalues)
%    =   0.2000    0.3020    0.4000
%        0.5020    0.6000    0.7020
%        0.8000    0.6000    0.2000
%        0.2000    0.2000    0.9020
% 
% 
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myRGBvalues = hex2rgb(myhexvalues,256)
%    =   51    77   102
%       128   153   179
%       204   153    51
%        51    51   230
% 
% HexValsAsACharacterArray = {'#334D66';'#8099B3';'#CC9933';'#3333E6'}; 
% rgbvals = hex2rgb(HexValsAsACharacterArray)
% 
% * * * * * * * * * * * * * * * * * * * * 
% Chad A. Greene, April 2014
%
% Updated August 2014: Functionality remains exactly the same, but it's a
% little more efficient and more robust. Thanks to Stephen Cobeldick for
% the improvement tips. In this update, the documentation now shows that
% the range may be set to 256. This is more intuitive than the previous
% style, which scaled values from 0 to 255 with range set to 255.  Now you
% can enter 256 or 255 for the range, and the answer will be the same--rgb
% values scaled from 0 to 255. Function now also accepts character arrays
% as input. 
% 
% * * * * * * * * * * * * * * * * * * * * 
% See also rgb2hex, dec2hex, hex2num, and ColorSpec. 
% 

%% Input checks:

assert(nargin>0&nargin<3,'hex2rgb function must have one or two inputs.') 

if nargin==2
    assert(isscalar(range)==1,'Range must be a scalar, either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end

%% Tweak inputs if necessary: 

if iscell(hex)
    assert(isvector(hex)==1,'Unexpected dimensions of input hex values.')
    
    % In case cell array elements are separated by a comma instead of a
    % semicolon, reshape hex:
    if isrow(hex)
        hex = hex'; 
    end
    
    % If input is cell, convert to matrix: 
    hex = cell2mat(hex);
end

if strcmpi(hex(1,1),'#')
    hex(:,1) = [];
end

if nargin == 1
    range = 1; 
end

%% Convert from hex to rgb: 

switch range
    case 1
        rgb = reshape(sscanf(hex.','%2x'),3,[]).'/255;

    case {255,256}
        rgb = reshape(sscanf(hex.','%2x'),3,[]).';
    
    otherwise
        error('Range must be either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end

end