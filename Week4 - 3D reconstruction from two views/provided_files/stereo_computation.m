function [ disparity_map ] = stereo_computation(left_image, right_image, minimum_disparity, maximum_disparity, window_size, matching_cost)

    % Check if image size is the same for both left and right images
    if size(left_image, 1) == size(right_image, 1) ...
            && size(left_image, 2) == size(right_image, 2)

        % Compute windows depending on window_size (check if even or odd)
        if (mod(window_size,2) == 0)
            wl = window_size/2; 
            wr = window_size/2-1;
        else
            wl = floor(window_size/2);
            wr = floor(window_size/2);
        end
        
        % Slide a window along the same line in the right image and compare
        % its content to that of the reference window in the left image.
        
        best_cost = [];
        disparity_map = [];
        
        for r = wl+1:size(left_image,1)-wr-1 %if img h=307 => 11:297
            for c = wl+1:size(left_image,2)-wr-1 %if img w=461 => 11:451
                
                w_im_1 = left_image(r-wl:r+wr, c-wl:c+wr); % window on left image
                
                for curr_disp = minimum_disparity:maximum_disparity
                    if(c+curr_disp-wl >= 1 && c+(wr+curr_disp) <= size(left_image,2))
                        w_im_2 = right_image(r-wl:r+wr, c+curr_disp-wl:c+curr_disp+wr); % window on right image
                        
                        % Compute ssd cost
                        if strcmp(matching_cost, 'ssd')
                            
                            if isempty(best_cost) && isempty(disparity_map)
                                disp('-- Computing ssd cost.')
                                best_cost = 99999.0 * ones(size(left_image));
                                disparity_map = zeros(size(left_image));
                            end
                            
                            cost = mean2(abs(w_im_1 - w_im_2).^2);
                        
                            % Look for the lowest
                            if cost < best_cost(r,c)
                                best_cost(r,c) = cost;
                                disparity_map(r,c) = curr_disp;
                            end
                        
                        % Compute ncc cost
                        elseif strcmp(matching_cost, 'ncc')
                            
                            if isempty(best_cost) && isempty(disparity_map)
                                disp('-- Computing ncc cost.')
                                best_cost = -1 * ones(size(left_image));
                                disparity_map = zeros(size(left_image));
                            end
   
                            cost = abs(mean2((w_im_1 - mean2(w_im_1)).*(w_im_2 - mean2(w_im_2))) / (std2(w_im_1) * std2(w_im_2)));
                            
                            % Look for the highest
                            if cost > best_cost(r,c)
                                best_cost(r,c) = cost;
                                disparity_map(r,c) = curr_disp;
                            end
                        
                        % Compute biweight cost
                        elseif strcmp(matching_cost, 'bilateral_weights')
                            
                            if isempty(best_cost) && isempty(disparity_map)
                                disp('-- Computing bilateral weight cost.')
                                best_cost = 99999.0 * ones(size(left_image));
                                disparity_map = zeros(size(left_image));
                            end
                            
                            cost = bilateral_weights_cost(w_im_1, w_im_2);
                            
                            % Look for the highest
                            if cost < best_cost(r,c)
                                best_cost(r,c) = cost;
                                disparity_map(r,c) = curr_disp;
                            end
                            
                        end
                        
                    end
                end
            end
        end
        
    else
        % Return an error value
        disparity_map = -99999.0;
        disp('Error: image size do not match.')
    end
end