   X = 1;
   Y = DR_CNO_CT.amp;
   bin_width = 1;
   spread_width = 0.1;



    val_min = min(Y);
    val_max = max(Y);
    Y = sort(Y);
    spread_rows = round((val_max - val_min)/bin_width);
    bin_idx = 0;
    swarm_x_temp = zeros(numel(Y),spread_rows);
    swarm_y_temp = zeros(numel(Y),spread_rows);
    
    for bi = val_min : bin_width: val_min+5*bin_width% : val_max
        bin_idx = bin_idx + 1;
        num_bin = 0;
        current_y = zeros(1,numel(Y));
        
        for yi = 1:numel(Y)
            if Y(yi) >= bi && Y(yi) < bi+bin_width
                num_bin = num_bin + 1;
                current_y(num_bin) = Y(yi);
            end
        end
        %if more than 1 values exist in this bin, spread.
         
        current_y = nonzeros(current_y);
        current_y_temp = zeros(1,num_bin);
        current_x = zeros(1,num_bin);    
        if num_bin > 1
            center_x = (num_bin+1)/2;
            center_y = current_y(end);
            current_y_temp(round(center_x)) = center_y;

            for xi = 1:num_bin
                if xi <= center_x
                    current_x(xi) = X(1)-(center_x - xi)*spread_width;
                else
                    current_x(xi) = X(1)+(xi - center_x)*spread_width;
                end
                
                swarm_x_temp(1:numel(current_x),bin_idx) = current_x;
            end
            %sort y as in gaussian distribution
            for yii = 1:2:(num_bin-1)
                 current_y_temp((yii+1)/2) = current_y(yii);
                 if num_bin > 2
                     current_y_temp(num_bin-(yii+1)/2+1) = current_y(yii+1);
                 end
            end
            
            swarm_y_temp(1:numel(current_y_temp),bin_idx) = current_y_temp;    
            
        elseif num_bin == 1
            swarm_x_temp(1,bin_idx) = X(1);
            swarm_y_temp(1,bin_idx) = current_y(1);
        end
    end
    
    swarm_X = nonzeros(swarm_x_temp);
    swarm_Y = nonzeros(swarm_y_temp);