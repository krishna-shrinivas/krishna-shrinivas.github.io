function generate_threshold_plot(FISH_data,IF_data,threshold_IF,normalize_flag,pos)
pos_x= pos(1);
pos_y = pos(2);
thresh_min = threshold_IF(1);
thresh_max = threshold_IF(2);
z_stacks = size(FISH_data{1},3);
count =0;  mean_IF = [];
xy_drift= 3;
subplot(2,2,1);
avg_intensity = zeros(size(FISH_data,2),z_stacks);
for i =1:1:size(FISH_data,2)
    
    for zpos =1:1:z_stacks
        for l=xy_drift*-1:1:xy_drift
            for m=xy_drift*-1:1:xy_drift
                avg_intensity(i,zpos) = avg_intensity(i,zpos)+double(IF_data{i}(pos_x+l,pos_y+m,zpos));
            end
        end
    end
    avg_intensity(i,:) = avg_intensity(i,:)/(2*xy_drift+1)^2;
    p = min(double(IF_data{i}(pos_x,pos_y,1:z_stacks)));
    if( (max(double(reshape(IF_data{i}(pos_x,pos_y,:),z_stacks,1))/p(1)) > thresh_min) && (max(double(reshape(IF_data{i}(pos_x,pos_y,:),z_stacks,1))/p(1)) < thresh_max))
        if normalize_flag
            plot(double(reshape(IF_data{i}(pos_x,pos_y,:),z_stacks,1))/p(1)); hold on; count= count+1;
            mean_IF(count,:) = double(reshape(IF_data{i}(pos_x,pos_y,:),z_stacks,1))/p(1);
        else
            plot(avg_intensity(i,:)); hold on; count= count+1;
            mean_IF(count,:) = double(avg_intensity(i,:));
        end
            
    end
end

title('Foci-resolved IF signal');
xlabel('Z coordinate');
ylabel('IF intensity');
subplot(2,2,2);
plot(mean(mean_IF)); hold on;
title('<IF> signal');
xlabel('Z coordinate');
ylabel('IF intensity');
disp(['No of foci with fold-change int > threshold is ' num2str(count-1)]);
avg_intensity = zeros(size(FISH_data,2),z_stacks);

subplot(2,2,3);
count =0;
mean_FISH = [];

for i =1:1:size(FISH_data,2)
    
    for zpos =1:1:z_stacks
        for l=xy_drift*-1:1:xy_drift
            for m=xy_drift*-1:1:xy_drift
                avg_intensity(i,zpos) = avg_intensity(i,zpos)+double(FISH_data{i}(pos_x+l,pos_y+m,zpos));
            end
        end
    end
    avg_intensity(i,:) = avg_intensity(i,:)/(2*xy_drift+1)^2;
    
    p = min(double(FISH_data{i}(pos_x,pos_y,1:z_stacks)));
    p_IF = min(double(IF_data{i}(pos_x,pos_y,1:z_stacks)));

    if( (max(double(reshape(IF_data{i}(pos_x,pos_y,:),z_stacks,1))/p_IF(1)) > thresh_min) && (max(double(reshape(IF_data{i}(pos_x,pos_y,:),z_stacks,1))/p_IF(1)) < thresh_max))
        if normalize_flag
            semilogy(double(reshape(FISH_data{i}(pos_x,pos_y,:),z_stacks,1))/p(1)); hold on; count= count+1;
            mean_FISH(count,:) = double(reshape(FISH_data{i}(pos_x,pos_y,:),z_stacks,1))/p(1);
        else
            plot(avg_intensity(i,:)); hold on; count= count+1;
            mean_FISH(count,:) = double(avg_intensity(i,:));
        end
    end
end
title('Foci-resolved FISH signal');
xlabel('Z coordinate');
ylabel('FISH intensity');

subplot(2,2,4);
plot(mean(mean_FISH)); hold on;
title('<FISH> signal');
xlabel('Z coordinate');
ylabel('Average intensity');
set(gcf, 'Position', get(0, 'Screensize'));


end