function [Xproj, Pproj] = factorization_method(x)
    % See lecture 8, page 50/51

    % 2. Normalize the set of points in each image
    for idx = 1:size(x, 2);
        [x_norm{idx}, H{idx}] = normalise2dpts(x{idx});
    end
    
    % 3. Initialize all projective depths (lambda)
    lambda = ones(length(x), size(x{1}, 2));
    for idx = 1:size(x, 2)
        F{idx} = fundamental_matrix(x{idx}, x{1});
        [U, D, V] = svd(F{idx});
        e{idx} = V(:,3) / V(3,3);
    end
    for i = 1:size(x, 2)
        for j = 1:size(x{1},2)
            num = x{1}(:, j)'*F{i}*cross(e{i}, x{i}(:,j));
            denom = norm(cross(e{i}, x{i}(:,j))).^2*lambda(1, j);
            lambda(i,j) = num/denom;
        end
    end

    % Initialize d and d_old
    % being d_old the distance in the previous iteration.
    d_old = 9e+10; d = 9e-10;
    
    % Initialize x_measurement, lambda multiplies this value to
    % obtain the measurement matrix M
    x_meas = [];
    
    for i=1:size(x, 2)
        x_meas = [x_meas; x_norm{i}];
    end
    
    % use the given convergence criterion
    while double(abs(d - d_old)/d) > 0.1

        %  4. Alternate rescaling the rows of the depth matrix ? to
        % have unit norm and the columns of ? to have unit norm until ?
        % stops changing significantly (usually two loops).
        for it = 1:2
            for i = 1:size(lambda,1)
                lambda(i,:) = lambda(i,:)/norm(lambda(i,:));
            end
            for i = 1:size(lambda,2)
                lambda(:,i) = lambda(:,i)/norm(lambda(:,i));
            end
        end

        prev = [];
        for i = 1:size(x, 2);
            prev = [prev; lambda(i, :); lambda(i, :); lambda(i, :)];
        end
        
        lambda = prev;
        d_old = d;

        % 5. Build the measurement matrix M
        M = lambda .* x_meas;
        
        % 6. Determine the SVD of M
        [U, D, V] = svd(M);

        % 7. Let PM = UD4 and XM = V_4'
        Pproj = U * D(:, 1:4);
        Xproj = V(:, 1:4)';

        % 8. If sum_i sum_j d() converges then stop, otherwise go to step 4.
        d = sum(sum((x_meas - Pproj * Xproj).^2));
        
        P = {};
        P{1} = Pproj(1:3,:);
        P{2} = Pproj(4:6,:);
        lambda = [];
        for i=1:length(P);
            aux = P{i} * Xproj;
            lambda = [lambda; aux(3,:)];
        end
    end
    % The last step is to unnormalize the camera matrices
    
    Pproj = [];
    for i = 1:size(H, 2)
        Pproj = [Pproj; mldivide(H{i}, P{i})];
    end

end