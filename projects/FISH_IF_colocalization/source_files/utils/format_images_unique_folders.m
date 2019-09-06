% function format_images_unique_folders(uniq)

directory_name = '/media/krishna/VERBATIM/RNA_Fish_Analysis/3D_IRF_analysis/images/new/';
common_name = 'V65_9_diff_Med1_FR_Nanog_G';
len_c = length(common_name);
p=dir(directory_name);
i=3;
while ( i <= length(p))
    
    local_files  = p(i).name;
    ind = strfind(local_files,common_name);
    if ~isempty(ind)
        uniq_folder_name = local_files(1+len_c+1:1+len_c+2);
        sub_folder = [directory_name 'Ediff' uniq_folder_name];
        if ~exist(sub_folder, 'dir') && ~isempty(sub_folder)
            mkdir(sub_folder);
        end
        file_name_complete = [directory_name local_files];
        status=movefile(file_name_complete,sub_folder);
        
    end
    
    i=i+1;
    
end
