function [ F ] = fundamental_matrix( x1, x2 )
% fundamental_matrix. Function in charge of computing 
% fundamental matrix from a set of 8 pairs of points
    
% Normalizate points
    [x1, H1] = normalise2dpts(x1);
    [x2, H2] = normalise2dpts(x2);
    
% Homogeneous coordinates
    x1(1,:) = x1(1,:)./x1(3,:);
    x1(2,:) = x1(2,:)./x1(3,:);
    x1(3,:) = 1;
    
    x2(1,:) = x2(1,:)./x2(3,:);
    x2(2,:) = x2(2,:)./x2(3,:);
    x2(3,:) = 1;

% u and v coordinates 
    u1 = x1(1,:)./x1(3,:);
    v1 = x1(2,:)./x1(3,:);

    u2 = x2(1,:)./x2(3,:);
    v2 = x2(2,:)./x2(3,:);
    
% Compute W matrix, such that Wf = 0
% Create matrix W from x1 and x2 correspondences
    W = [];
    for i = 1:size(x1,2)
       tmp = [u1(i)*u2(i) v1(i)*u2(i) u2(i) u1(i)*v2(i) v1(i)*v2(i) v2(i) u1(i) v1(i) 1];
       W = [W; tmp];        
    end
    
% Compute the SVD of matrix W=UDV'
[~,~,V] = svd(W);

% Compose fundamental matrix F_rank3
F_rank3 = reshape(V(:,end),3,3)';

% Compute the SVD of fundamental matrix F_rank3
[U,D,V] = svd(F_rank3);

% Remove last singular value of D to create D_rank2
D_rank2 = D;
D_rank2(end,end) = 0;

% Re-compute matrix F = U D_rank2 V' (rank 2)
F = U * D_rank2 * V';

% De-normalize F 
F = H2' * F * H1;
end

