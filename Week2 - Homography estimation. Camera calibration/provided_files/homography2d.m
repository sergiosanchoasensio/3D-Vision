function H = homography2d(x,xp)
% Rebust homography computation by means of the 
% Direct Linear Transformation (DLT), by group 7.

% Step 1: normalize.
% This method is not invariant to similarity transformations of the image.
% Thus some coordinate systems are in some way better than others for
% computing a 2D homography. 
% Solution to this problem: to apply a method of normalization of the
% data (consisting of translation and scaling of image coordinates) before
% applying the DLT algorithm.
[x,t1] = Normalise2DPts(x);
[xp,t2] = Normalise2DPts(xp);


% Step 2: Ah = 0. The vector h is in the null space of A.
xp_x = xp(1,:);
xp_y = xp(2,:);
xp_w = xp(3,:);
A = [];
% Computing homographies: from Lecture 2, page 42/60
for i=1:size(x,2)
    A = [A; zeros(3,1)'     -xp_w(i)*x(:,i)'   xp_y(i)*x(:,i)'; ...
            xp_w(i)*x(:,i)'   zeros(3,1)'     -xp_x(i)*x(:,i)'];
end

% Step 3: compute the SVD of A. The unit singular vector corresponding to the
% smallest singular value is the solution h. A = UDV' Then h is the last column of V.
[~,~,v] = svd(A);
h = reshape(v(:,9),3,3)';

% Step 4: Desnormalization
H = inv(t2) * h * t1;
end


% Normalized DLT algorithm: from Lecture 2, page 54/60
function [p2, T] = Normalise2DPts(p)
    % This normalizing transformation will diminish the effect of the arbitrary
    % selection of origin and scale in the coordinate frame of the image, and
    % will mean that the combined algorithm is invariant to a similarity
    % transformation of the image.
    
    % For x: Compute a similarity transformation T ,
    % consisting of a translation and scaling, that takes points x_i to a new
    % set of points ~x_i such that the centroid of the points ~x_i is the
    % coordinate origin, and their average distance from the origin is sqrt(2).
    
    % For x': Normalization of x' : Compute a similar transformation T' for the
    % points in the second image, transforming points x_i' to a new set of
    % points ~x_i'
    
    p(1,:) = p(1,:)./p(3,:);
    p(2,:) = p(2,:)./p(3,:);
    p(3,:) = 1;
    
    c = mean(p(1:2,:)')';
    p2(1,:) = p(1,:)-c(1);
    p2(2,:) = p(2,:)-c(2);
    
    scale = sqrt(2)/mean(sqrt(p2(1,:).^2 + p2(2,:).^2));
    
    T = [scale   0   -scale*c(1)
         0     scale -scale*c(2)
         0       0      1      ];
    
    p2 = T*p;
end
