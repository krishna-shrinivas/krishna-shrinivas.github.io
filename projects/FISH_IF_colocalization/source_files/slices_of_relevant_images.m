% Function that reads in FISH and IF files to generate output centered on
% on various FISH foci and saves all the captured foci, as well as the IF
% information in all the stacks for any given FISH foci
%
%   The various parameters are
%
%       directory_name ==   directory under which the files are stored under
%                           with each replicate under a different subfolder
%
%       flag == 1 if passed FISH area & intensity threshold is used
%              0 if default threshold of 8 pixel^2, 2.1 multiplier is used
%
%       threshold_multiplier == passed intensity multiplier threshold
%
%       FISH_threshold == passed FISH area threshold
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
%

function slices_of_relevant_images(directory_name,output_dir_name,input_params)


%Lists all sub-directories in the images folder which
% are characterise as Expt_no_1/Image_1, Expt_no/Image_2
p=dir(directory_name);
subdir = {};

FISH_data = {};
IF_data = {};
COM_data = [];

%Pixel parameters
xpixel = input_params.xpixel;
zpixel = input_params.zpixel;

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
%probe csv folder if exists:

if ((input_params.curate_called_foci) || (~input_params.automatically_call_foci))
   p = dir(input_params.csv_folder);
   count = 1;
   for id=1:1:length(uniq_folder_id)
       i=3;
       flag =0;
       while ((i<=length(p)) && (flag==0))
            indx = strfind(p(i).name,uniq_folder_id{id});
            if ~isempty(indx)
                curated_foci_name{id} = [input_params.csv_folder p(i).name];
                flag =1;
            end
            i=i+1;
       end
    end
end


%   Thresholds for max and min areas of FISH foci per stack
FISH_area_threshold = 1;
FISH_max_area_threshold = 200;
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


if input_params.flag
    FISH_area_threshold = input_params.FISH_threshold;
end


foci_ac= 1;
image_data_set = {};

for file =1:1:length(IF_channel)
    n_images=numel(imfinfo(IF_channel{file}));
    count =1;
    Image_store = {};
    Fish_store= [];
    stats_relevant = [];
    filename = IF_channel{file};
    filename_FISH = FISH_channel{file};
    
    for stack = 1:1:n_images
        [X,map] = imread(filename,stack);
        X_store{file,stack} =X;
    end
    for stack=1:1:n_images
        % Stack of BRD4
        
        [XF,mapF] = imread(filename_FISH,stack);
        XF_store{file,stack} = XF;
        
        if ~input_params.flag
            mean_intensity = mean(mean(XF));
            
            bwF= imbinarize(XF,2.1*mean_intensity/65535);
        else
            mean_intensity = mean(mean(XF));
            bwF= imbinarize(XF,input_params.threshold_multiplier*mean_intensity/65535);
        end
        
        bwF = bwareaopen(bwF,FISH_area_threshold);
        store_data{file,stack} = bwF;
        [BF,LF] = bwboundaries(bwF,'noholes');
        
        
        stats = regionprops(LF,'BoundingBox','Area','Centroid', 'MajorAxisLength','MinorAxisLength','PixelIdxList');
        %
        
        B_rel = [];
        stat_rel = [];
        size_box = input_params.size_box;
        if ~isempty(stats)
            
            for b_rel = 1:length(BF)
                if (stats(b_rel).Area >= FISH_area_threshold) && (stats(b_rel).Area <FISH_max_area_threshold)
                    Centroid = stats(b_rel).Centroid;
                    
                    if (floor(Centroid(1)) - size_box >1) && (floor(Centroid(1)) + size_box<size(XF,1)) && (floor(Centroid(2)) + size_box<size(XF,1)) && (floor(Centroid(2)) - size_box>1)
                        flag_end_stack_points = 0;
                        if (stack==1)
                            flag_end_stack_points =1;
                        elseif (stack<= n_images-max_offset)
                            for local_stack = 1:1:size(offset_range,2)
                                Image_store{count}(:,:,local_stack) = X_store{file,stack+offset_range(local_stack)}(floor(Centroid(2)) - size_box:floor(Centroid(2)) + size_box,floor(Centroid(1)) - size_box:floor(Centroid(1)) + size_box);
                            end
                        else
                            flag_end_stack_points =1;
                        end
                        
                        if ~flag_end_stack_points
                            Fish_store(:,:,count) =  XF_store{file,stack}(floor(Centroid(2)) - size_box:floor(Centroid(2)) + size_box,floor(Centroid(1)) - size_box:floor(Centroid(1)) + size_box);
                            stats_relevant(count,1:2) = Centroid;
                            stats_relevant(count,3) = stack;
                            stats_relevant(count,4) = stats(b_rel).Area;
                            stats_relevant(count,5) = mean(XF_store{file,stack}(stats(b_rel).PixelIdxList));
                            count = count+1;
                        end
                    end
                end
            end
        end
    end
    
    
    if ~isempty(Image_store)
        for offset = 1:1:size(offset_range,2)
            
            IF_store = [];
%             curated_foci_name{file} = [gene_folder '/Foci_calls/' pair_of_genes '_' p(file+2).name(end-1:end) '.' FISH_name '.csv'];
            
            
            for i =1:1:size(Fish_store,3)
                IF_store(:,:,i) = Image_store{i}(:,:,offset);
            end
            
            if input_params.automatically_call_foci
                [COM]= find_foci_centroid(Fish_store,IF_store,stats_relevant,input_params);
                if input_params.curate_called_foci
                    [COM] = curate_foci_COM(curated_foci_name{file},COM,input_params);
                end
            else
                called_foci = csvread(curated_foci_name{file},1,4);
                called_foci(:,4:6) = [];
                COM = [];
                COM(:,1:2) = ceil(called_foci(:,1:2)/xpixel);
                COM(:,3) = round(called_foci(:,3)/zpixel);
            end
            
            size_box_xy = input_params.size_box;
            size_box_z = floor((size_box_xy*xpixel)/zpixel);
            for foci = 1:1:size(COM,1)
                Centroid = COM(foci,:);
                stack_centroid = Centroid(3);
                
                
                
                if (stack_centroid - size_box_z >1) && (stack_centroid +size_box_z < n_images - offset_range(offset))
                    if (floor(Centroid(1)) - size_box >1) && (floor(Centroid(1)) + size_box<size(XF,1)) && (floor(Centroid(2)) + size_box<size(XF,1)) && (floor(Centroid(2)) - size_box>1)
                        
                        for stack = 1:1:2*size_box_z +1
                            FISH_data{foci_ac}(:,:,stack) = XF_store{file,round(stack_centroid) -size_box_z -1+stack+offset_range(offset)}(floor(Centroid(2)) - size_box:floor(Centroid(2)) + size_box,floor(Centroid(1)) - size_box:floor(Centroid(1)) + size_box);
                            IF_data{foci_ac}(:,:,stack) = X_store{file,round(stack_centroid) -size_box_z -1+stack+offset_range(offset)}(floor(Centroid(2)) - size_box:floor(Centroid(2)) + size_box,floor(Centroid(1)) - size_box:floor(Centroid(1)) + size_box);
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
end


for offset=1:1:size(offset_range,2)
    raw_data = [parentFolder '/Combined_data/' 'All_data_' num2str(offset_range(offset)) '.mat'];
    combined_IRF_figure = [parentFolder '/Combined_data/' 'All_IRF_' num2str(offset_range(offset))];
    combined_average_figure = [parentFolder '/Combined_data/' 'All_average_' num2str(offset_range(offset))];

    [temp_parent,temp] =fileparts(raw_data);
    if ~exist(temp_parent, 'dir') && ~isempty(temp_parent)
        mkdir(temp_parent);
    end
    if ~isempty(FISH_data)
        save([temp_parent '/Total_data.mat'],'FISH_IRF','IF_IRF','FISH_data','IF_data','COM_data','image_data_set','input_params');
        generate_IRF_collapse_plot(FISH_IRF,IF_IRF,input_params.ind_flag,combined_IRF_figure);
        generate_average_images(IF_data,FISH_data,combined_average_figure,IF_name,FISH_name);
    end
end
end
%

