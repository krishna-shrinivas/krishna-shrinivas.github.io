path_to_real_data = 'John_Imaging/18_05_04_Trim28RNAFISH/trim28_oct4_size_box_30/Combined_data/'
path_to_random_data = 'John_Imaging/18_05_04_Trim28RNAFISH/trim28_oct4_size_box_30_random/Combined_data/'

% path_to_real_data = '/media/krishna/VERBATIM/RNA_FISH_Analysis/3D_IRF_Analysis/Pol II studies/20180505_Med1_Pol2_RNAFISH/trim28_polIICTDS2P_med1/trim28_polIICTDS2P_auto-calledxy10z10_size_box_30/Combined_data/'
IF_name  = 'Oct4';
FISH_name = 'Trim28';

real_data_loc = [path_to_real_data 'Total_data.mat'];
random_data_loc = [path_to_random_data 'Total_data.mat'];

FISH_color = [1 0 1];
IF_color = [0 1 0];
RD = load(real_data_loc);
Random_data =load(random_data_loc);
image_offset_loc =1;

name = [path_to_real_data  FISH_name '_' IF_name '_random_control_8_10'];
production_average_images(RD.IF_data,RD.FISH_data,Random_data.IF_data,name,IF_name,FISH_name,FISH_color,IF_color,[8,10])

% figure_name = [path_to_real_data  FISH_name '_' IF_name '_new_colorscheme_production_IRF_image'];
% production_quality_IRF(RD.FISH_IRF,RD.IF_IRF,0,figure_name,FISH_name,IF_name,FISH_color,IF_color)
