function [] = auto_call_random_foci(directory_name,random_foci_dir,input_params)
%   Function calls randomly located centroid position within nuclei,
%   identified by the DAPI channel.
%   Input parameters to the function are:
%       o directory_name - link to /.../images
%       o random_foci_dir - Output directory where called foci are written,
%       typically '/...../Random_foci_auto'
%       o input_params - Input parameters - refer params.m for more info.

p=dir(directory_name);
subdir = {};

%   Pixel parameters
xpixel = input_params.xpixel;
zpixel = input_params.zpixel;

%   Stores DNA FISH parameters
for i=3:1:length(p)
    subdir{i-2} = [directory_name p(i).name '/'];
    if isdir(subdir{i-2})
        uniq_folder_id{i-2} = p(i).name(end-1:end);
        local_files  = dir(subdir{i-2});
        for j=3:1:length(local_files)
            if strfind(local_files(j).name,'405')
                DNA_channel{i-2} =[subdir{i-2} local_files(j).name];
                output_DNA_file{i-2} = [ random_foci_dir 'auto_called_random_foci_' uniq_folder_id{i-2}];
            end
        end
        cell_classifier_DNA(DNA_channel{i-2},output_DNA_file{i-2},input_params,40);
        
    end
end
end