function [IRF] = convert_3D_image_IRF(Data,input_params)
%   Converts 3D matrix of image data into a 1D radial function about centroid
center_point = [ceil(size(Data,1)/2) ceil(size(Data,2)/2) ceil(size(Data,3)/2)];

xsize = size(Data,1);
ysize = size(Data,2);
zsize = size(Data,3);
count =1;

xpixel  = input_params.xpixel;
zpixel = input_params.zpixel;

for i =1:1:xsize
    for j=1:1:ysize
        for k =1:1:zsize
            
            dist_center(count)= sqrt(xpixel*xpixel*((i-center_point(1))^2+(j-center_point(2))^2) + zpixel*zpixel*(k-center_point(3))^2);
            Intensity(count) = double(Data(i,j,k));
            count = count+1;
        end
    end
end


[sort_dist_center,index_sort]=sort(dist_center);
[sort_Intensity] = Intensity(index_sort);


   
d= 1;
while(d <= length(sort_dist_center))
    repeats = find(sort_dist_center== sort_dist_center(d));
    if length(repeats)~=1
        sort_dist_center(repeats(2:end)) = [];
        sort_Intensity(d) = mean(sort_Intensity([d, repeats]));        
        sort_Intensity(repeats(2:end)) = [];
    end
    d=d+1;
end

IRF = [];
IRF.dist = sort_dist_center;
IRF.intensity = double(sort_Intensity);








end