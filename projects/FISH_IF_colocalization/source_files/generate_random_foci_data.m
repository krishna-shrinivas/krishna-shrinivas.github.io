%   Function that reads in FISH, IF, and random nuclear spot files to generate output centered on
%   on various *random* FISH foci and saves all the captured foci, as well as the IF
%   information in all the stacks for any given FISH foci
%
%   The various parameters are
%
%       directory_name ==   directory under which the files are stored under
%                           with each replicate under a different subfolder
%
%       output_dir_name ==  directory to which output files are written
%                           In particular, the there will be a separate
%                           subfolder for each replicate and one labeled
%                           'Combined_data' with all data together. Within
%                           each replicate folder, there will be many
%                           sub-folders:
%                               Data - Foci calls and MATLAB data
%                               Images_offset_x -  Average images z-offset by
%                                                  x slices b/w FISH & IF
%                               IRF_offset_x -     IRF of average images
%                                                  with _av_0 for <IRF> and
%                                                  _av_1 IRF of <image>
%
%       random_nuclear_foci ==  Folder with random nuclear foci

function generate_random_foci_data(directory_name,output_dir_name,random_nuclear_foci,input_params)


%Lists all sub-directories in the images folder which
% are characterise as Expt_no_1/Image_1, Expt_no/Image_2
p=dir(directory_name);
subdir = {};
COM_data = [];
%   Stores IF and FISH information different cells
%   for each separate replicate
count = 1;
for i=3:1:length(p)
    
    
    subdir{i-2} = [directory_name p(i).name '/'];
    if (isdir(subdir{i-2}) && isempty(strfind(subdir{i-2},'.')))
        
        local_files  = dir(subdir{i-2});
        for j=3:1:length(local_files)
            if strfind(local_files(j).name,input_params.IF)
                IF_channel{count} =[subdir{i-2} local_files(j).name];
                uniq_folder_id{count} = p(i).name(end-1:end);
            elseif strfind(local_files(j).name,input_params.FISH)
                FISH_channel{count} =[subdir{i-2} local_files(j).name];

            end 
        end
        count= count+1;
    end
end

%   Data structure that keeps track of all the foci across experiments
overall_images = [];

%   Range of z-stack offsets from min_offset to max_offset
max_offset = 0;
min_offset = 0;
offset_range = min_offset:1:max_offset;

[outputdir filen]= fileparts(directory_name);
[parentFolder filen] =fileparts(outputdir);
gene_folder = parentFolder;


[temp pair_of_genes] = fileparts(parentFolder);

FISH_name = input_params.FISH_name;
IF_name = input_params.IF_name;

parentFolder = [parentFolder '/' output_dir_name];
if ~exist(parentFolder, 'dir') && ~isempty(parentFolder)
    mkdir(parentFolder);
end
overall_images.I = cell(max_offset-min_offset+1,1);
overall_images.S = cell(max_offset-min_offset+1,1);
overall_images.F = cell(max_offset-min_offset+1,1);

image_data_set = {};

foci_ac= 1;


xyumtopixelconversion = 1/input_params.xpixel;
zumtopixelconversion = 1/input_params.zpixel;

xpixel = input_params.xpixel;
zpixel = input_params.zpixel;


%% Store random foci data
p = dir([gene_folder '/' random_nuclear_foci]);
count = 1;
for id=1:1:length(IF_channel)
    i=3;
    flag =0;
    while ((i<=length(p)) && (flag==0))
        indx = strfind(p(i).name,uniq_folder_id{id});
        indx2 = strfind(p(i).name,'txt');
        % to ensure that file has a txt in it
        if (~isempty(indx) && ~isempty(indx2))
            curated_foci_name{id} = [gene_folder '/' random_nuclear_foci '/' p(i).name];
            flag =1;
        end
        i=i+1;
    end
end



for file =1:1:length(IF_channel)
    
    randomly_called_foci{file} = csvread(curated_foci_name{file},1,1);
    randomly_called_foci{file}(:,1:2) = randomly_called_foci{file}(:,1:2) * xyumtopixelconversion;
    randomly_called_foci{file}(:,3) = randomly_called_foci{file}(:,3) * zumtopixelconversion;
    randomly_called_foci{file} = floor(randomly_called_foci{file});
    
    n_images=numel(imfinfo(IF_channel{file}));
    count =1;
    Image_store = {};
    Fish_store= [];
    stats_relevant = [];
    filename = IF_channel{file};
    filename_FISH = FISH_channel{file};
    
    % Reading all the IF and FISH images
    for stack = 1:1:n_images
        [X,map] = imread(filename,stack);
        X_store{file,stack} =X;
        
        [XF,mapF] = imread(filename_FISH,stack);
        XF_store{file,stack} = XF;
        
    end
    
    
    for offset = 1:1:size(offset_range,2)
        
        
        COM = [];
        COM(:,1:2) = randomly_called_foci{file}(:,1:2);
        COM(:,3) = randomly_called_foci{file}(:,3);
        
        size_box = input_params.size_box;
        size_box_z = floor((size_box*xpixel)/zpixel);
        for foci = 1:1:size(COM,1)
            Centroid = COM(foci,:);
            stack_centroid = Centroid(3);
            
            
            if (stack_centroid - size_box_z >1) && (stack_centroid +size_box_z < n_images - offset_range(offset))
                if (floor(Centroid(1)) - size_box >1) && (floor(Centroid(1)) + size_box<size(XF,1)) && (floor(Centroid(2)) + size_box<size(XF,1)) && (floor(Centroid(2)) - size_box>1)
                    
                    for stack = 1:1:2*size_box_z +1
                        FISH_data{foci_ac}(:,:,stack) = XF_store{file,stack_centroid -size_box_z -1+stack+offset_range(offset)}(floor(Centroid(2)) - size_box:floor(Centroid(2)) + size_box,floor(Centroid(1)) - size_box:floor(Centroid(1)) + size_box);
                        IF_data{foci_ac}(:,:,stack) = X_store{file,stack_centroid -size_box_z -1+stack+offset_range(offset)}(floor(Centroid(2)) - size_box:floor(Centroid(2)) + size_box,floor(Centroid(1)) - size_box:floor(Centroid(1)) + size_box);
                    end
                    [FISH_IRF{foci_ac}] = convert_3D_image_IRF(FISH_data{foci_ac},input_params);
                    [IF_IRF{foci_ac}] = convert_3D_image_IRF(IF_data{foci_ac},input_params);
                    image_data_set{foci_ac} = filename;
                    COM_data(foci_ac,:) = COM(foci,:);
                    foci_ac =foci_ac +1;
                    
                end
            end
            
        end
        
    end
    
    
    
    
    
    
end


for offset=1:1:size(offset_range,2)
    raw_data = [parentFolder '/Combined_data/' 'All_data_' num2str(offset_range(offset)) '.mat'];
    combined_IRF_figure = [parentFolder '/Combined_data/' 'All_IRF_' num2str(offset_range(offset))];
    combined_average_figure = [parentFolder '/Combined_data/' 'All_average_' num2str(offset_range(offset))];
    
    [temp_parent,temp] =fileparts(raw_data);
    if ~exist(temp_parent, 'dir') && ~isempty(temp_parent)
        mkdir(temp_parent);
    end
    save([temp_parent '/Total_data.mat'],'FISH_IRF','IF_IRF','FISH_data','IF_data','COM_data','image_data_set','input_params');
    generate_IRF_collapse_plot(FISH_IRF,IF_IRF,input_params.ind_flag,combined_IRF_figure)
    generate_average_images(IF_data,FISH_data,combined_average_figure,IF_name,FISH_name)
end




end
%

