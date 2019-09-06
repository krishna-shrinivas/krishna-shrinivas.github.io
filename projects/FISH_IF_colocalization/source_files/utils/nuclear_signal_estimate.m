
function[IF_mean] = nuclear_signal_estimate(filename,IF_file)
close all;


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
IF_signal = [];
rand_store = [];
count = 1;
for stack=1:1:length(stack_to_check)
    bw = im2bw(XDNA_store{file,stack_to_check(stack)},(mean(p{stack_to_check(stack)}) + threshold*std(p{stack_to_check(stack)}))/65535);
    bw = bwareaopen(bw,100);
    [B,L] = bwboundaries(bw,'noholes');
    stats = regionprops(L,'Area','Centroid', 'MajorAxisLength','MinorAxisLength','PixelIdxList');
    k_rel = [];

    [X_IF,map] = imread(IF_file,stack_to_check(stack));
    
    for k = 1:length(B)
        boundary = B{k};
        if length(boundary) > 200
           IF_signal = [IF_signal; X_IF(stats(k).PixelIdxList)];    
        end
    end



end

IF_mean = mean(IF_signal);

end
