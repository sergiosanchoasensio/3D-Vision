%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 2: Image mosaics

clear all; clc; close all;

addpath('sift');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Compute image correspondences

%% Open images
% imargb = imread('Data/llanes/llanes_a.jpg');
% imbrgb = imread('Data/llanes/llanes_b.jpg');
% imcrgb = imread('Data/llanes/llanes_c.jpg');

imargb = imread('Data/castle_int/0016_s.png');
imbrgb = imread('Data/castle_int/0015_s.png');
imcrgb = imread('Data/castle_int/0014_s.png');

%imargb = imread('Data/aerial/site13/frame00000.png');
%imbrgb = imread('Data/aerial/site13/frame00002.png');
%imcrgb = imread('Data/aerial/site13/frame00003.png');

ima = sum(double(imargb), 3) / 3 / 255;
imb = sum(double(imbrgb), 3) / 3 / 255;
imc = sum(double(imcrgb), 3) / 3 / 255;

% imargb = double(imread('Data/aerial/site22/frame_00001.tif'));
% imbrgb = double(imread('Data/aerial/site22/frame_00018.tif'));
% imcrgb = double(imread('Data/aerial/site22/frame_00030.tif'));
% ima = imargb;
% imb = imbrgb;
% imc = imcrgb;

%% Compute SIFT keypoints
[points_a, desc_a] = sift(ima, 'Threshold', 0.01);
[points_b, desc_b] = sift(imb, 'Threshold', 0.01);
[points_c, desc_c] = sift(imc, 'Threshold', 0.01);

figure;
imshow(imargb);%image(imargb)
hold on;
plot(points_a(1,:), points_a(2,:),'+y');
figure;
imshow(imbrgb);%image(imbrgb);
hold on;
plot(points_b(1,:), points_b(2,:),'+y');
figure;
imshow(imcrgb);%image(imcrgb);
hold on;
plot(points_c(1,:), points_c(2,:),'+y');

%% Match SIFT keypoints 

% between a and b
matches_ab = siftmatch(desc_a, desc_b);
figure;
plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), matches_ab, 'Stacking', 'v');

% between b and c
matches_bc = siftmatch(desc_b, desc_c);
figure;
plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), matches_bc, 'Stacking', 'v');

% between c and b
matches_cb = siftmatch(desc_c, desc_b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Compute the homography (DLT algorithm) between image pairs

%% Compute homography (normalized DLT) between a and b, play with the homography
th = 3;
xab_a = [points_a(1:2, matches_ab(1,:)); ones(1, length(matches_ab))];
xab_b = [points_b(1:2, matches_ab(2,:)); ones(1, length(matches_ab))];
[Hab, inliers_ab] = ransac_homography_adaptive_loop(xab_a, xab_b, th, 1000); % ToDo: complete this function

figure;
plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), ...
    matches_ab(:,inliers_ab), 'Stacking', 'v');

vgg_gui_H(imargb, imbrgb, Hab);


%% Compute homography (normalized DLT) between b and c, play with the homography
xcb_c = [points_c(1:2, matches_cb(1,:)); ones(1, length(matches_cb))];
xcb_b = [points_b(1:2, matches_cb(2,:)); ones(1, length(matches_cb))];
[Hcb, inliers_cb] = ransac_homography_adaptive_loop(xcb_c, xcb_b, th, 1000);

figure;
plotmatches(imc, imb, points_c(1:2,:), points_b(1:2,:), ...
    matches_cb(:,inliers_cb), 'Stacking', 'v');

vgg_gui_H(imbrgb, imcrgb, Hcb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Build the mosaic

corners = [-400 1200 -100 650];
I = [1 0 0; 0 1 0; 0 0 1];

% b is the image in the middle
iwb = apply_H_v2(imbrgb, I, corners);   % ToDo: complete the call to the function
% a to b
iwa = apply_H_v2(imargb, Hab, corners);    % ToDo: complete the call to the function
% c to b
iwc = apply_H_v2(imcrgb, Hcb, corners);    % ToDo: complete the call to the function

figure;
imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));axis off;
title('Mosaic A-B-C by normalized DLT');

% ToDo: compute the mosaic with castle_int images
% ToDo: compute the mosaic with aerial images set 13
% ToDo: compute the mosaic with aerial images set 22
% ToDo: comment the results in every of the four cases: say why it works or
%       does not work

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Refine the homography with the Gold Standard algorithm
% Homography ab
x = xab_a(1:2, inliers_ab);  %ToDo: set the non-homogeneous point coordinates of the 
xp = xab_b(1:2, inliers_ab); % point correspondences we will refine with the geometric method

Xobs = [ x(:) ; xp(:) ];     % The column vector of observed values (x and x')
P0 = [ Hab(:) ; x(:) ];      % The parameters or independent variables

Y_initial = gs_errfunction( P0, Xobs ); % ToDo: create this function that we need to pass to the lsqnonlin function
% NOTE: gs_errfunction should return E(X) and not the sum-of-squares E=sum(E(X).^2)) that we want to minimize. 
% (E(X) is summed and squared implicitly in the lsqnonlin algorithm.) 

err_initial = sum( sum( Y_initial.^2 ));
options = optimset('Algorithm', 'levenberg-marquardt');
P = lsqnonlin(@(t) gs_errfunction(t, Xobs), P0, [], [], options);
Hab_r = reshape( P(1:9), 3, 3 );
f = gs_errfunction( P, Xobs ); % lsqnonlin does not return f
err_final = sum( sum( f.^2 ));
% we show the geometric error before and after the refinement
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initial, err_final);
%% See differences in the keypoint locations
% ToDo: compute the points xhat and xhatp which are the correspondences
% returned by the refinement with the Gold Standard algorithm
xhat = P(10:end);
xhat = reshape(xhat,2,[]);
one = ones(1,size(xhat,2));
xhathomog = [xhat; one];

xhatp = Hab_r*xhathomog;
xhatp = xhatp ./ repmat(xhatp(end,:), size(xhatp,1), 1);

figure;
imshow(imargb);%image(imargb);
hold on;
plot(x(1,:), x(2,:),'+y');
plot(xhat(1,:), xhat(2,:),'+c');
figure;
imshow(imbrgb);%image(imbrgb);
hold on;
plot(xp(1,:), xp(2,:),'+y');
plot(xhatp(1,:), xhatp(2,:),'+c');
%%  Homography bc
% ToDo: refine the homography bc with the Gold Standard algorithm
x = xcb_c(1:2, inliers_cb);  %ToDo: set the non-homogeneous point coordinates of the 
xp = xcb_b(1:2, inliers_cb); % point correspondences we will refine with the geometric method

Xobs = [ x(:) ; xp(:) ];     % The column vector of observed values (x and x')
P0 = [ Hcb(:) ; x(:) ];      % The parameters or independent variables

Y_initial = gs_errfunction( P0, Xobs ); % ToDo: create this function that we need to pass to the lsqnonlin function
% NOTE: gs_errfunction should return E(X) and not the sum-of-squares E=sum(E(X).^2)) that we want to minimize. 
% (E(X) is summed and squared implicitly in the lsqnonlin algorithm.) 

err_initial = sum( sum( Y_initial.^2 ));
options = optimset('Algorithm', 'levenberg-marquardt');
P = lsqnonlin(@(t) gs_errfunction(t, Xobs), P0, [], [], options);
Hcb_r = reshape( P(1:9), 3, 3 );
f = gs_errfunction( P, Xobs ); % lsqnonlin does not return f
err_final = sum( sum( f.^2 ));
% we show the geometric error before and after the refinement
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initial, err_final);
%% See differences in the keypoint locations
% ToDo: compute the points xhat and xhatp which are the correspondences
% returned by the refinement with the Gold Standard algorithm
xhat = P(10:end);
xhat = reshape(xhat,2,[]);
one = ones(1,size(xhat,2));
xhathomog = [xhat; one];

xhatp = Hcb_r*xhathomog;
xhatp = xhatp ./ repmat(xhatp(end,:), size(xhatp,1), 1);

figure;
imshow(imbrgb);%image(imbrgb);
hold on;
plot(x(1,:), x(2,:),'+y');
plot(xhat(1,:), xhat(2,:),'+c');
figure;
imshow(imcrgb);%image(imcrgb);
hold on;
plot(xp(1,:), xp(2,:),'+y');
plot(xhatp(1,:), xhatp(2,:),'+c');
%% Build mosaic
corners = [-400 1200 -100 650];
iwb = apply_H_v2(imbrgb, I, corners); % ToDo: complete the call to the function
iwa = apply_H_v2(imargb, Hab_r, corners); % ToDo: complete the call to the function
iwc = apply_H_v2(imcrgb, Hcb_r, corners); % ToDo: complete the call to the function
figure;
imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));axis off;
title('Mosaic A-B-C by Gold Standard algorithm');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Calibration with a planar pattern

clear all; close all;

%% Read template and images.
T     = imread('Data/calib/template.jpg');
I{1}  = imread('Data/calib/graffiti1.tif');
I{2}  = imread('Data/calib/graffiti2.tif');
I{3}  = imread('Data/calib/graffiti3.tif');
%I{4}  = imread('Data/calib/graffiti4.tif');
%I{5}  = imread('Data/calib/graffiti5.tif');
Tg = sum(double(T), 3) / 3 / 255;
Ig{1} = sum(double(I{1}), 3) / 3 / 255;
Ig{2} = sum(double(I{2}), 3) / 3 / 255;
Ig{3} = sum(double(I{3}), 3) / 3 / 255;
N = length(I);

%% Compute keypoints.
fprintf('Computing sift points in template... ');
[pointsT, descrT] = sift(Tg, 'Threshold', 0.05);
fprintf(' done\n');
points = cell(N,1);
descr = cell(N,1);
for i = 1:N
    fprintf('Computing sift points in image %d... ', i);
    [points{i}, descr{i}] = sift(Ig{i}, 'Threshold', 0.05);
    fprintf(' done\n');
end
%% Match and compute homographies.
H = cell(N,1);
for i = 1:N
    % Match against template descriptors.
    fprintf('Matching image %d... ', i);
    matches = siftmatch(descrT, descr{i});
    fprintf('done\n');
    % Fit homography and remove outliers.
    x1 = pointsT(1:2, matches(1, :));
    x2 = points{i}(1:2, matches(2, :));
    
    % Convert to homog. coords
    x1_h = [x1; ones(1, size(x1,2))];
    x2_h = [x2; ones(1, size(x2,2))];

    H{i} = 0;
    [H{i}, inliers] =  ransac_homography_adaptive_loop(x1_h, x2_h, 3, 1000);
    % Plot inliers.
    figure;
    plotmatches(Tg, Ig{i}, pointsT(1:2,:), points{i}(1:2,:), matches(:, inliers));
    % Play with the homography
    vgg_gui_H(T, I{i}, H{i});
end

%% Compute the Image of the Absolute Conic
% Image of the Absolute Conic, which is usually written as w. 
% For computational pruposes, it will be simpler to first find w
% and then recover the camera matrix K from it.
% We know that hi'*w*hj = Vij'*x. vij'*x = 0 and (vii'-vjj')*x = 0

v = [];
for n=1:N

    H_i = H{n};
    
    i = 1; 
    j = 2;
    v12_T = [H_i(1,i)*H_i(1,j), H_i(1,i)*H_i(2,j) + H_i(2,i)*H_i(1,j), ...
        H_i(1,i)*H_i(3,j) + H_i(3,i)*H_i(1,j), H_i(2,i)*H_i(2,j), ...
        H_i(2,i)*H_i(3,j) + H_i(3,i)*H_i(2,j), H_i(3,i)*H_i(3,j)];
            
    i = 1; 
    j = 1;
    v11_T = [H_i(1,i)*H_i(1,j), H_i(1,i)*H_i(2,j) + H_i(2,i)*H_i(1,j), ...
        H_i(1,i)*H_i(3,j) + H_i(3,i)*H_i(1,j), H_i(2,i)*H_i(2,j), ...
        H_i(2,i)*H_i(3,j) + H_i(3,i)*H_i(2,j), H_i(3,i)*H_i(3,j)];
        
    i = 2; 
    j = 2;
    v22_T = [H_i(1,i)*H_i(1,j), H_i(1,i)*H_i(2,j) + H_i(2,i)*H_i(1,j), ...
        H_i(1,i)*H_i(3,j) + H_i(3,i)*H_i(1,j), H_i(2,i)*H_i(2,j), ...
        H_i(2,i)*H_i(3,j) + H_i(3,i)*H_i(2,j), H_i(3,i)*H_i(3,j)];
    
    v1 = v12_T;
    v2 = v11_T - v22_T;
    
    v = [v;v1;v2];  
end

[~, ~, V] = svd(v);

omega = V(:,end);

w = [omega(1), omega(2), omega(3); 
     omega(2), omega(4), omega(5); 
     omega(3), omega(5), omega(6)]; % ToDo

%% Recover the camera calibration.
 K = chol(inv(w),'upper');
 
%%
        % ToDo: in the report make some comments related to the obtained internal
        %       camera parameters and also comment their relation to the image size
        
        % Add comments 

%% Compute camera position and orientation.
R = cell(N,1);
t = cell(N,1);
P = cell(N,1);
figure;hold;
for i = 1:N
    
    H_i = H{i};
    
    % ToDo: compute r1, r2, and t{i}
    r1 = inv(K) * H_i(:,1) / norm(inv(K) * H_i(:,1));
    r2 = inv(K) * H_i(:,2) / norm(inv(K) * H_i(:,2));
    t{i} = inv(K) * H_i(:,3) / norm(inv(K) * H_i(:,1));
    
    % Solve the scale ambiguity by forcing r1 and r2 to be unit vectors.
    s = sqrt(norm(r1) * norm(r2)) * sign(t{i}(3));
    r1 = r1 / s;
    r2 = r2 / s;
    t{i} = t{i} / s;
    R{i} = [r1, r2, cross(r1,r2)];
    
    % Ensure R is a rotation matrix
    [U S V] = svd(R{i});
    R{i} = U * eye(3) * V';
   
    P{i} = K * [R{i} t{i}];
    plot_camera(P{i}, 800, 600, 200);
end
%%
        % ToDo: in the report explain how the optical center is computed in the
        %       provided code

        % Add comments 

%%
[ny,nx] = size(T);
p1 = [0 0 0]';
p2 = [nx 0 0]';
p3 = [nx ny 0]';
p4 = [0 ny 0]';
% Draw planar pattern
vgg_scatter_plot([p1 p2 p3 p4 p1], 'g');
% Paint image texture
surface('XData',[0 nx; 0 nx],'YData',[0 0; 0 0],'ZData',[0 0; -ny -ny],'CData',T,'FaceColor','texturemap');
colormap(gray);
axis equal;

%% Plot a static camera with moving calibration pattern.
figure; hold;
plot_camera(K * eye(3,4), 800, 600, 200);

% ToDo: complete the call to the following function with the proper
%       coordinates of the image corners in the new reference system
corners = [p1 p2 p3 p4 p1];
corners(end+1, : ) = 1;
for i = 1:N  
    calibrationPoints = inv(K) * (P{i} * corners);
    vgg_scatter_plot(calibrationPoints, 'r');
end

%% Augmented reality: Plot some 3D points on every camera.
[Th, Tw] = size(Tg);
cube = [0 0 0; 1 0 0; 1 0 0; 1 1 0; 1 1 0; 0 1 0; 0 1 0; 0 0 0; 0 0 1; 1 0 1; 1 0 1; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 0 1; 0 0 0; 1 0 0; 1 0 0; 1 0 1; 1 0 1; 0 0 1; 0 0 1; 0 0 0; 0 1 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 1 0; 0 0 0; 0 1 0; 0 1 0; 0 1 1; 0 1 1; 0 0 1; 0 0 1; 0 0 0; 1 0 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 1 0 1; 1 0 1; 1 0 0 ]';
X = (cube - .5) * Tw / 4 + repmat([Tw / 2; Th / 2; -Tw / 8], 1, length(cube));
X(end+1,:) = 1;
for i = 1:N
    figure; colormap(gray);
    imagesc(Ig{i});
    hold on;
    x = euclid(P{i} * X);
    vgg_scatter_plot(x, 'g');
end

% ToDo: change the virtual object, use another 3D simple geometric object like a pyramid
[Th, Tw] = size(Tg);
x = [ 0 1 1 1 1 0 0 0 sqrt(2)/2 1 sqrt(2)/2 0 sqrt(2)/2 1 sqrt(2)/2];
y = [ 0 0 0 0 0 0 0 0 -1 0 -1 0 -1 0 -1]; 
z = [ 0 0 0 1 1 1 1 0 sqrt(2)/2 0 sqrt(2)/2 1 sqrt(2)/2 1 sqrt(2)/2];

pyramid = [];
for i = 1:size(x,2)
    pyramid = [pyramid; [x(i) y(i) z(i)]];
end

theta = degtorad(-90);
Rx_theta = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
%Ry_theta = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
%Rz_theta = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

pyramid = pyramid*Rx_theta;

pyramid = pyramid';
X = (pyramid - repmat([.5; 0; .5], 1, size(pyramid,2))) * Tw / 4 + ...
    repmat([Tw / 2; Th / 2; -Tw / 8], 1, size(pyramid,2));
X(end+1,:) = 1;
for i = 1:N
    figure; colormap(gray);
    imagesc(Ig{i});
    hold on;
    x = euclid(P{i} * X);
    vgg_scatter_plot(x, 'r');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. OPTIONAL: Add a logo to an image using the DLT algorithm
clear all; close all; clc;
imargb = imread('Data/logo/timesSquare.jpg');
imbrgb = imread('Data/logo/bw.jpg');

ima = sum(double(imargb), 3) / 3 / 255;
imb = sum(double(imbrgb), 3) / 3 / 255;

% Points image 1
p1 = [1488 384 1]';
p2 = [1957 267 1]';
p3 = [1972 832 1]';
p4 = [1500 844 1]';
x1 = [p1 p2 p3 p4];

% Points image 2
[ny,nx] = size(imb);
p1 = [0 0 1]';
p2 = [nx 0 1]';
p3 = [nx ny 1]';
p4 = [0 ny 1]';
x2 = [p1 p2 p3 p4];

H = homography2d(x2, x1);

[ny,nx] = size(ima);
corners = [0 nx 0 ny];
I = [1 0 0; 0 1 0; 0 0 1];

% a is the image in the middle
iwa = apply_H_v2(imargb, I, corners);   % ToDo: complete the call to the function
% b to a
iwb = apply_H_v2(imbrgb, H, corners);    % ToDo: complete the call to the function

e = find(iwb);
iwa(e) = 0;

figure;
imshow(max(iwb, iwa));
title('Added image into another image by normalized DLT');
