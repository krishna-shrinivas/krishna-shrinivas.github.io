function [COM] = curate_foci_COM(curated_foci_csv,COM,input_params)
%   Function that curates auto-called foci vs manually called foci
xpixel = input_params.xpixel;
zpixel = input_params.zpixel;
COM(:,1:2) = COM(:,1:2)*xpixel;
COM(:,3) = COM(:,3)*zpixel;

dist_threshold =input_params.distance_threshold;
called_foci = csvread(curated_foci_csv,1,4);
called_foci(:,4:6) = [];
dist_called_curated = [];
remove_foci  = [];
for foci =1:1:size(COM,1)
    for curated_foci = 1:1:size(called_foci,1)
        dist_called_curated(foci,curated_foci) = sqrt(sum((COM(foci,1:3) - called_foci(curated_foci,:)).^2));
    end
    if min(dist_called_curated(foci,:)) > dist_threshold
        remove_foci = [remove_foci foci];
    end
end
if ~isempty(remove_foci)
    COM(remove_foci,:) = [];

end


COM(:,1:2) = ceil(COM(:,1:2)/xpixel);
COM(:,3) = ceil(COM(:,3)/zpixel);

end