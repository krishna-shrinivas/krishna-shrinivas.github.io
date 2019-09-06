%   Parametrizes input with typical image-analysis parameters
%   Certain parameters and flag variables can be set within this
%   function. The list is:
%       o threshold_multiplier - Multiplicative threshold for intensity
%       thresholding
%       o FISH_threshold - Minimum area (in pixels) for thresholding
%       o xpixel - Pixel length (um) along x/y directions
%       o zpixel - Pixel length (um) along z stack
%       o distance_threshold - Pixel displacement (um) across z-stacks of
%       centroid that is permissible for stitching
%       o automatically_call_foci - Flag variable that decides whether (1)
%       or not (0) automatic foci detection should happen. Default setting
%       is ON (1).
%       o curate_called_foci - Flag variable that decides whether (1)
%       or not (0) automatically called foci are curated across a manually
%       provided csv list of foci. Default setting is OFF (0).
%       o volume_threshold - Minimum volume (um^3) for accepting foci.
%       Default setting is 0.05 um^3 <-- experience for FISH data in mESCs.
%       o size_box - Half length of box around FISH centroid (in pixels)
%       for which the IF data is gathered. Default is 25 pixels ~ 1.5 um.
%       o random_auto_call - Flag variable that decides whether (1)
%       or not (0) to call random foci based on DAPI channel data. Default
%       setting is ON (1) - requires 405 channel image data.
%
%
%   This function completely parametrizes the KEY variable - input_params,
%   which is used to peform the image-colocalization analysis. The
%   function does not return any data - instead generating all the key
%   output files under the output_dir_name.

function params(directory_name,output_dir_name,in_params)

input_params = [];

%   Flag dictates which internal threshold to employ
%   Use threshold specified below
input_params.flag = 1;

if isempty(in_params)
    input_params.FISH = '561';
    input_params.IF = '488';

    input_params.FISH_name = 'mir290';
    input_params.IF_name = 'med1';
else
    input_params.FISH = in_params.FISH_channel;
    input_params.IF = in_params.IF_channel;

    input_params.FISH_name = in_params.FISH_name;
    input_params.IF_name = in_params.IF_name;
end

%   Passed minimum intensity multiplier and area (pixel^2)
%   threshold for FISH

input_params.threshold_multiplier = 2.0;
input_params.FISH_threshold = 12.0;

%   um pixel resolution for microscope

input_params.xpixel = 0.0572;
input_params.zpixel = 0.2;

%   This is the distance threshold for centroid stitching across multiple
%   z-stacks ( in um )
input_params.distance_threshold = 0.75;

%   Setting automatically_call_foci =1 allows the code to call FISH foci, and
%   setting it to 0 turns off automated foci calling
input_params.automatically_call_foci = 1;

%   Setting curate_called_foci=1 curates the automatically called foci with
%   the manually called foci stored in files stored under
%   'C:\.....\mir290_brd4\$CSV_FOLDER', where csv_folder is the folder with
%   called foci files that is defined below. Setting curate_called_foci=0
%   does not curate the called data.

input_params.curate_called_foci = 0;

%   Volumetric threshold for accepting a called FISH foci, in units of um^3
input_params.volume_threshold = 0.05;

%   csv_folder points to the folder at 'C:\.....\mir290_brd4\csv_folder\'
%   where the files for manually called foci are stored.
input_params.csv_folder = 'Foci_calls';
input_params.csv_folder = [fileparts(fileparts(directory_name)) '/' input_params.csv_folder '/'];

%   This is the half the length of the cube of image data that is stored,
%   centered on the centroid of the FISH foci. In units of xy pixels. The
%   size of the z plane is calculated as (size_box)* (x_pixel)/(z_pixel).
input_params.size_box = 25;

input_params.ind_flag = 0;


%   This is the name of the folder under which the randomly called foci
%   centroids are stored.'C:\.....\mir290_brd4\$random_nuclear_foci'.
random_nuclear_foci = 'Random_foci';

%   Auto call random foci- suite still under development
input_params.random_auto_call = 1;

output_dir_name = [ output_dir_name '_size_box_' num2str(input_params.size_box) ];

% File-name for Nuclear-IF signal - debugging ON for next two lines
% nuclear_out_file = [fileparts(fileparts(directory_name)) '/' output_dir_name '/nuclear_IF_average.csv'];
% normalized_data_collection(directory_name,nuclear_out_file,input_params);

%   This calls the functions where the majority of analysis calculations
%   take place.
slices_of_relevant_images(directory_name,output_dir_name,input_params);

% If the random foci are called and stored in the folder, then the function
% will call the random foci analysis.
random_exist = isdir([fileparts(fileparts(directory_name)) '/' random_nuclear_foci '/']);
if random_exist
    output_random_dir_name = [ output_dir_name '_random'];
    generate_random_foci_data(directory_name,output_random_dir_name,random_nuclear_foci,input_params)

elseif input_params.random_auto_call
    random_foci_dir = [fileparts(fileparts(directory_name)) '/' random_nuclear_foci '_auto/'];
    if ~isdir(random_foci_dir)
        mkdir([fileparts(fileparts(directory_name)) '/' random_nuclear_foci '_auto/']);
        auto_call_random_foci(directory_name,random_foci_dir,input_params);
    end
    output_random_dir_name = [ output_dir_name '_random'];
    generate_random_foci_data(directory_name,output_random_dir_name,[random_nuclear_foci '_auto/'],input_params)

end
%% In this section, once the raw data is generated and stored, production quality images
%% are generated and saved in the output_directory with RED for FISH and green for IF.
real_data_loc = [fileparts(fileparts(directory_name)) '/' output_dir_name  '/Combined_data/Total_data.mat'];
%   Saves important foci statistics
save_foci_data_post_analysis(real_data_loc,input_params);

FISH_color = [1 0 0];
IF_color = [0 1 0];
RD = load(real_data_loc);
z_stack_center = ceil(size(RD.FISH_data{1},3)/2);
name = [fileparts(real_data_loc) '/'  input_params.FISH_name '_' input_params.IF_name '_production_average_image'];
figure_name = [fileparts(real_data_loc) '/'   input_params.FISH_name '_' input_params.IF_name '_production_IRF_image'];

if ((~random_exist) && (~input_params.random_auto_call))
    production_average_images(RD.IF_data,RD.FISH_data,[],name,input_params.IF_name,input_params.FISH_name,FISH_color,IF_color,[z_stack_center-5,z_stack_center+5])
    production_quality_IRF(RD.FISH_IRF,RD.IF_IRF,[],0,figure_name,input_params.FISH_name,input_params.IF_name,FISH_color,IF_color)

else
    random_data_loc = [fileparts(fileparts(directory_name)) '/' output_random_dir_name  '/Combined_data/Total_data.mat'];
    save_foci_data_post_analysis(random_data_loc,input_params);

    Random_store =load(random_data_loc);
    production_average_images(RD.IF_data,RD.FISH_data,Random_store.IF_data,name,input_params.IF_name,input_params.FISH_name,FISH_color,IF_color,[z_stack_center-5,z_stack_center+5])
    production_quality_IRF(RD.FISH_IRF,RD.IF_IRF,Random_store.IF_IRF,0,figure_name,input_params.FISH_name,input_params.IF_name,FISH_color,IF_color)


end


end
