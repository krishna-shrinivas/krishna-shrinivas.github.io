function[] = cell_classifier_DNA(filename,output_file,input_params,N_rand)
%   Function that takes in a DAPI-channel to identify nuclei.
%   Typically called from auto_call_random_foci.m Key parameters are:
%       o filename - link to /.../images/E01/whatever_image_05_405.TIF
%       o output_file - Output directory where called foci are written for
%       above file typically '/...../Random_foci_auto/auto_called_random_foci_05.csv'
%       o input_params - Input parameters - refer params.m for more info.
%       o N_rand - Number of foci to call per image-set
%
%
%   General DAPI segmentation is based on an intensity threshold and
%   average tracking of stack intensity that has been empirically
%   determined to work well for DAPI stains in mESCs

close all;

%   Set draw=1 to one to debug individual DAPI files
draw = 0;
file = 1;
n_images=numel(imfinfo(filename));
p = {};
for stack = 1:1:n_images
    [X,map] = imread(filename,stack);
    XDNA_store{file,stack} =X;
    p{stack}=double(reshape(X,length(X)^2,1));
    mean_data(stack) = mean(p{stack});
end
cell_id_thre = 0.1;
stack_to_check = find(mean_data> 0.1*(max(mean_data)-min(mean_data))+ min(mean_data));
threshold =0.25;
list_of_points = [];
rand_store = [];
count = 1;
for stack=1:1:length(stack_to_check)
    bw = im2bw(XDNA_store{file,stack_to_check(stack)},(mean(p{stack_to_check(stack)}) + threshold*std(p{stack_to_check(stack)}))/65535);
    bw = bwareaopen(bw,100);
    [B,L] = bwboundaries(bw,'noholes');
    stats = regionprops(L,'Area','Centroid', 'MajorAxisLength','MinorAxisLength','PixelIdxList');
    k_rel = [];
    for k = 1:length(B)
        boundary = B{k};
        if length(boundary) > 200
            list_of_points = [ list_of_points ;stats(k).PixelIdxList];
            k_rel = [k_rel k];
        end
    end

    

    random_points = list_of_points(randperm(length(list_of_points),5));
    [I,J] = ind2sub([length(bw), length(bw)],random_points);
    if draw
        figure;
        imshow(XDNA_store{file,stack_to_check(stack)},[]);
        hold on
        for k = 1:length(k_rel)
            boundary = B{k_rel(k)};
            if length(boundary) > 200
                plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2); hold on;
                scatter(J,I,'blue','filled','d'); hold on;
            end
        end
        drawnow;
    end
    rand_store(count:count+4,1:4) = [[count:count+4]',J,I, stack_to_check(stack)*ones(5,1);];
    count = count+5;
end
rand_store(:,2:3) =rand_store(:,2:3)*0.0572;
rand_store(:,4) =rand_store(:,4)*0.2;
subset_pick = randperm(length(rand_store(:,1)),N_rand);
subset_return = [];
subset_return(:,2:4) = rand_store(subset_pick,2:4);
subset_return(:,1) = 1:1:N_rand';
T = array2table(subset_return);
T.Properties.VariableNames = {'Foci_number','X_centroid','Y_centroid','Z_centroid'};
writetable(T, output_file);

end
