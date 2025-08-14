function [ ] = forward_particle_start_accessibilities(init_r_list, impacted, grouping, q_over_m)

impacted = ~impacted;

y_vals = unique(init_r_list(:,2));
z_vals = unique(init_r_list(:,3));
Z = zeros(length(y_vals), length(z_vals));

for i = 1:length(y_vals)
    for j = 1:length(z_vals)
        count = 0;
        total = 0;
        for k = 1:size(init_r_list, 1)
            if init_r_list(k, 2) == y_vals(i) && init_r_list(k, 3) == z_vals(j)
                if impacted(k)
                    count = count + 1;
                end
                total = total + 1;
            end
        end
        if total > 0
            Z(i, j) = count / total * 100;  % now: row = lat, col = long
        else
            Z(i, j) = NaN;
        end
    end
end

if grouping
    Z_new = zeros(length(y_vals)./grouping(1), length(z_vals)./grouping(2));
    for i = 1:length(y_vals)
        for j = 1:length(z_vals)
            inds = [ceil(i./grouping(1)), ceil(j./grouping(2))];
            Z_new(inds(1), inds(2)) = Z_new(inds(1), inds(2)) + Z(i, j);
        end
    end
    Z_new = Z_new ./ (grouping(1) * grouping(2));

    yinds = arrayfun(@(i) mean(y_vals(i:min(i+grouping(1)-1,end))), ...
                     1:grouping(1):length(y_vals));
    zinds = arrayfun(@(j) mean(z_vals(j:min(j+grouping(2)-1,end))), ...
                     1:grouping(2):length(z_vals));

    figure;
    clf;
    imagesc(yinds, zinds, Z_new);    
    xlim([min(y_vals), max(y_vals)]);
    ylim([min(z_vals), max(z_vals)]);
    set(gca,'YDir','normal'); 
    xlabel('y (m)') 
    ylabel('z (m)')                                        
    colormap(jet);
    cb = colorbar;
    cb.Label.String = 'Impact Rates (%)';
    clim([0 100]);   
    title(sprintf('Starting Position Impact Rates, q/m = %d', q_over_m))
else
    figure;
    clf;
    imagesc(y_vals, z_vals, Z);              
    set(gca,'YDir','normal'); 
    xlabel('y (m)')
    ylabel('z (m)')                                         
    colormap(jet);
    cb = colorbar;
    cb.Label.String = 'Impact Rates (%)';
    clim([0 100]);
    title(sprintf('Starting Position Impact Rates, q/m = %d', q_over_m))
end


