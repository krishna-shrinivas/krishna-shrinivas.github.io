function []= save_foci_data_post_analysis(real_data_location,params)
%   Analysis function that saves foci_data post output data generation

RD = load(real_data_location);
[output_dir_name,b] =fileparts(real_data_location);

center = ceil(size(RD.FISH_data{1})/2);
xy_drift = 3;
z_drift = 1;
count = (2*xy_drift+1)*(2*xy_drift+1)*(2*z_drift+1);
average_intensity_roi = zeros(length(RD.FISH_data),2);
file_names ={};
for foci = 1:1:length(RD.FISH_data)
    
    for p=xy_drift*-1:1:xy_drift
        for q=xy_drift*-1:1:xy_drift
            
            for r =z_drift*-1:1:z_drift
                average_intensity_roi(foci,1) = average_intensity_roi(foci,1)+ double(RD.FISH_data{foci}(center(1)+p,center(2)+q,center(3)+r));
                average_intensity_roi(foci,2) = average_intensity_roi(foci,2) + double(RD.IF_data{foci}(center(1)+p,center(2)+q,center(3)+r));
                
            end
        end
    end
    file_names{foci} ='';
    [c,file_names{foci}]=fileparts(RD.image_data_set{foci});
end

average_intensity_roi = average_intensity_roi/count;
%         corrcoef(average_intensity_roi(:,2),average_intensity_roi(:,3))

output_data = [];
output_data(:,1:3) = RD.COM_data(:,1:3);
output_data(:,1:2) =output_data(:,1:2) * params.xpixel;
output_data(:,3) =output_data(:,3) * params.zpixel;

output_data(:,4:5) = average_intensity_roi;
T = array2table(output_data);
T(:,6) = array2table(file_names');
T.Properties.VariableNames ={'X_center','Y_center','Z_center','Average_FISH','Average_IF', 'Experimental_image_file'};

output_file = [output_dir_name '/' params.FISH_name '_' params.IF_name '_foci_statistics.csv'];
writetable(T,output_file);
close all;
%         generate_threshold_plot(RDM.FISH_data,RDM.IF_data,[0 26],0,[21 21]);
%         output_image_name = [FISH_name '_med1_all_foci.svg'];
%         saveas(gcf,output_image_name);
%         close all;

time_to_complete = toc;
end