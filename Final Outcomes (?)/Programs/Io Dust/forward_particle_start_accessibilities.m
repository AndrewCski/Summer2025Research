function [ ] = forward_particle_start_accessibilities(init_r_list, impacted, grouping, q_over_m, moon)

    % Function for plotting the accessibilities of particles started in a
    % uniform grid moving towards a moon. init_r_list is the position
    % vectors of the initial states of each particle, impacted is an array
    % indicating whether each particle impacted the moon or not, grouping
    % is a value indicating how particles should be binned together, if
    % at all, q_over_m is the q/m ratio of the particles, and moon is a
    % string indicating whether Europa or Ganymede is simulated.

impacted = ~impacted;

y_vals = unique(init_r_list(:,2));
z_vals = unique(init_r_list(:,3));
Z = zeros(length(y_vals), length(z_vals));

if moon == "Europa"
    moon_rad = 1560e3;
    moon_string = "(R_E)";
elseif moon == "Ganymede"
    moon_rad = 2631e3;
    moon_string = "(R_G)";
end

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
    imagesc(yinds./moon_rad, zinds./moon_rad, Z_new);    
    xlim([min(y_vals./moon_rad), max(y_vals./moon_rad)]);
    ylim([min(z_vals./moon_rad), max(z_vals./moon_rad)]);
    set(gca,'YDir','normal'); 
    xlabel(sprintf('y %s', moon_string)) 
    ylabel(sprintf('z %s', moon_string))                                        
    colormap(jet);
    cb = colorbar;
    cb.Label.String = 'Impact Rates (%)';
    clim([0 100]);   
    if all(q_over_m == q_over_m(1), 'all')
        title(sprintf('Starting Position Impact Rates, q/m = %d', q_over_m(1)))
    else
        title('Starting Position Impact Rates')
    end
else
    figure;
    clf;
    imagesc(y_vals./moon_rad, z_vals./moon_rad, Z);              
    set(gca,'YDir','normal'); 
    xlabel(sprintf('y %s', moon_string)) 
    ylabel(sprintf('z %s', moon_string))                                         
    colormap(jet);
    cb = colorbar;
    cb.Label.String = 'Impact Rates (%)';
    clim([0 100]);
    if all(q_over_m == q_over_m(1), 'all')
        title(sprintf('Starting Position Impact Rates, q/m = %d', q_over_m(1)))
    else
        title('Starting Position Impact Rates')
    end
end


