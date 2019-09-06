function write_foci_tiff_image(FISH_data,IF_data,output_file_IF,output_file_FISH)

%   Write individual foci into a long TIFF file
mean_FISH = FISH_data{1};
mean_IF = IF_data{1};
for i=2:1:size(FISH_data,2)
    mean_FISH = mean_FISH + FISH_data{i};
    mean_IF = mean_IF + IF_data{i};
end

mean_FISH = mean_FISH/size(FISH_data,2);
mean_IF = mean_IF/size(IF_data,2);

imwrite(mean_FISH(:,:,1),output_file_FISH);
imwrite(mean_IF(:,:,1),output_file_IF);

for j=2:1:size(FISH_data{i},3)
    imwrite(mean_FISH(:,:,j),output_file_FISH,'WriteMode','append');
    imwrite(mean_IF(:,:,j),output_file_IF,'WriteMode','append');
end


end