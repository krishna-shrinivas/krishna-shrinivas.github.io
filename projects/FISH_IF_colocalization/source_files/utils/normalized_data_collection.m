function [] = normalized_data_collection(directory_name,output_file,input_params)


p=dir(directory_name);
subdir = {};

%Pixel parameters
xpixel = input_params.xpixel;
zpixel = input_params.zpixel;

%   Stores DNA FISH + IF for nuclear counting parameters
for i=3:1:length(p)
    subdir{i-2} = [directory_name p(i).name '/'];
    if isdir(subdir{i-2})
        uniq_folder_id{i-2} = p(i).name(end-1:end);
        local_files  = dir(subdir{i-2});
        for j=3:1:length(local_files)
            if strfind(local_files(j).name,'405')
                DNA_channel{i-2} =[subdir{i-2} local_files(j).name];
            end
            if strfind(local_files(j).name,'488')
                IF_channel{i-2} = [subdir{i-2} local_files(j).name];
            end
        end

        IF_mean{i-2} = nuclear_signal_estimate(DNA_channel{i-2},IF_channel{i-2});
        
    end
end
all_data = [IF_channel' IF_mean'];
T = cell2table(all_data);
T.Properties.VariableNames = {'Image_file','Average_IF'};

if ~exist((fileparts(output_file)), 'dir') && ~isempty((fileparts(output_file)))
    mkdir((fileparts(output_file)));
end
writetable(T, output_file);

end