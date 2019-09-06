function [COM]= find_foci_centroid(Fish_store,IF_store,stats_relevant,input_params)
%   Function takes in FISH data and performs registration/thresholding across
%   the z-stacks to identify COM for each foci

N_Fish_store = Fish_store;
N_Image_store = IF_store;
N_stats = stats_relevant(:,1:5);


zthreshold = 10;
xythreshold = 10;

xpixel = input_params.xpixel;
zpixel = input_params.zpixel;

foci =1;
while (foci<= size(N_stats,1))
    COM = N_stats(foci,1:2);
    stack = N_stats(foci,3);
    store_to_remove = [];
    
    for j=foci+1:1:size(N_stats,1)
        dist = sqrt ( sum((COM(1,1:2)- N_stats(j,1:2)).^2));
        if (dist < xythreshold) && ( N_stats(j,3)-stack <=zthreshold)
            store_to_remove = [store_to_remove j];
        end
    end
    
    
    if (COM(1)>1000) && (COM(2)>1000)
        N_stats(foci,:) = [];
        N_Fish_store(:,:,foci) = [];
        N_Image_store(:,:,foci) = [];
        
    else
        N_stats(foci,1) = sum(N_stats([foci store_to_remove],1).*N_stats([foci store_to_remove],4))/sum(N_stats([foci store_to_remove],4));
        N_stats(foci,2) = sum(N_stats([foci store_to_remove],2).*N_stats([foci store_to_remove],4))/sum(N_stats([foci store_to_remove],4));
        N_stats(foci,3) = sum(N_stats([foci store_to_remove],3).*N_stats([foci store_to_remove],4))/sum(N_stats([foci store_to_remove],4));
        N_stats(foci,5) = sum(N_stats([foci store_to_remove],4) .* N_stats([foci store_to_remove],5))/sum(N_stats([foci store_to_remove],4));
        N_stats(foci,4) = sum(N_stats([foci store_to_remove],4));
        
        N_Fish_store(:,:,foci) = mean(N_Fish_store(:,:,[foci store_to_remove]),3);
        N_Image_store(:,:,foci) = mean(N_Image_store(:,:,[foci store_to_remove]),3);
        
        N_stats(store_to_remove,:) = [];
        N_Fish_store(:,:,store_to_remove) = [];
        N_Image_store(:,:,store_to_remove) = [];
        foci = foci+1;
        
    end
    
end



N_stats(:,6:7) = N_stats(:,1:2)*xpixel;
N_stats(:,8) = N_stats(:,3)*zpixel;
N_stats(:,9) =N_stats(:,4)*(xpixel*xpixel)*zpixel;

volume_threshold = input_params.volume_threshold;
foci =1;
while (foci<= size(N_stats,1))
    
    Vol_cur = N_stats(foci,9);
    if Vol_cur < volume_threshold
        N_stats(foci,:) = [];
        N_Fish_store(:,:,foci) = [];
        N_Image_store(:,:,foci) = [];
        foci= foci-1;
    end
    foci= foci+1;
end



COM = [];
COM(:,1:3) = N_stats(:,1:3);

%Volume of bodies
COM(:,4) = N_stats(:,9);
end