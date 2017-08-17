%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Lab 3: The geometry of two views 
% (optional: photo-sequencing)

addpath('sift'); % ToDo: change 'sift' to the correct path where you have the sift functions
clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Compute the fundamental matrix

% Two camera matrices for testing purposes
P1 = eye(3,4);
c = cosd(15); s = sind(15);
R = [c -s 0; s c 0; 0 0 1];
t = [.3 0.1 0.2]';
P2 = [R t];
n = 8;
X = [rand(3,n); ones(1,n)] + [zeros(2,n); 3 * ones(1,n); zeros(1,n)];
x1_test = P1 * X;
x2_test = P2 * X;

% Estimated fundamental matrix
% ToDo: create the following function that estimates F using the normalised 8 point algorithm
F_es = fundamental_matrix(x1_test, x2_test);

% Real fundamental matrix
% F=K2'^-1[t]RK1^-1. In our case, K1 = K2 = I beacuse P1 = I and P2 = [R t];
Tx = [0 -t(3) t(2); t(3) 0  -t(1); -t(2) t(1) 0];
F_gt = Tx * R; % ToDo: write the expression of the real fundamental matrix for P1 and P2

% Evaluation: these two matrices should be very similar
F_gt / norm(F_gt)
F_es / norm(F_es)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Robustly fit fundamental matrix

% Read images
im1rgb = imread('Data/0000_s.png');
im2rgb = imread('Data/0001_s.png');
im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;

% show images
figure;
subplot(1,2,1); imshow(im1rgb); axis image; title('Image 1');
subplot(1,2,2); imshow(im2rgb); axis image; title('Image 2');


%% Compute SIFT keypoints

% (make sure that the sift folder provided in lab2 is on the path)

[points_1, desc_1] = sift(im1, 'Threshold', 0.01);
[points_2, desc_2] = sift(im2, 'Threshold', 0.01);

%% Match SIFT keypoints between a and b
matches = siftmatch(desc_1, desc_2);
figure;
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches, 'Stacking', 'v');

% p1 and p2 contain the homogeneous coordinates of the matches
p1 = [points_1(1:2, matches(1,:)); ones(1, length(matches))];
p2 = [points_2(1:2, matches(2,:)); ones(1, length(matches))];

% ToDo: create this function (you can use as a basis 'ransac_homography_adaptive_loop.m')
[F, inliers] = ransac_fundamental_matrix(p1, p2, 2.0); 

% show inliers
figure;
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches(:,inliers), 'Stacking', 'v');
title('Inliers');

vgg_gui_F(im1rgb, im2rgb, F');


%% Plot some epipolar lines
number_of_inliers = 3;
points = inliers(randsample(size(inliers,2), number_of_inliers));

% choose three random indices
m1 = points(1);
m2 = points(2);
m3 = points(3);

ip1 = p1(:,points);  %pick coordinates of the inliers points
ip2 = p2(:,points);     
l1 = F'*ip2; % epipolar lines in image 2
l2 = F*ip1;% epipolar lines in image 1 % ToDo

% image 1 (plot the three points and their corresponding epipolar lines)
figure;
imshow(im1rgb);
hold on;
plot(p1(1, m1), p1(2, m1), '+g');
plot_homog_line(l1(:, 1)); %indices of the points m1,m2,m3

plot(p1(1, m2), p1(2, m2), '+g');
plot_homog_line(l1(:, 2)); %indices of the points m1,m2,m3

plot(p1(1, m3), p1(2, m3), '+g');
plot_homog_line(l1(:, 3)); %indices of the points m1,m2,m3

% image 2 (plot the three points and their corresponding epipolar lines)
figure;
imshow(im2rgb);
hold on;
plot(p2(1, m1), p2(2, m1), '+g');
plot_homog_line(l2(:, 1)); %indices of the points m1,m2,m3

plot(p2(1, m2), p2(2, m2), '+g');
plot_homog_line(l2(:, 2)); %indices of the points m1,m2,m3

plot(p2(1, m3), p2(2, m3), '+g');
plot_homog_line(l2(:, 3)); %indices of the points m1,m2,m3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. OPTIONAL: Photo-sequencing with aerial images

% In this part we will compute a simplified version of the algorithm
% explained in the Photo-sequencing paper. 
% Since we do not have two images
% taken from roughly the same viewpoint at two different time instants we
% will manually pick a dynamic point corresponding to a point in a van 
% (identified by index 'idx_car_I1') and the projection of its 3D trajectory 
% in the reference image. Then we will compute the projection (to the reference image) 
% of three points on this 3D trajectory at three different time instants 
% (corresponding to the time when the three other provided images where taken). 

clear all;

% Read images
im1rgb = imread('Data/frame_00000.tif');
im2rgb = imread('Data/frame_00001.tif');
im3rgb = imread('Data/frame_00002.tif');
im4rgb = imread('Data/frame_00003.tif');

im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;
im3 = sum(double(im3rgb), 3) / 3 / 255;
im4 = sum(double(im4rgb), 3) / 3 / 255;

% show images
figure;
subplot(2,2,1); imshow(im1rgb); axis image; title('Image 1');
subplot(2,2,2); imshow(im2rgb); axis image; title('Image 2');
subplot(2,2,3); imshow(im3rgb); axis image; title('Image 3');
subplot(2,2,4); imshow(im4rgb); axis image; title('Image 4');

% Compute SIFT keypoints
[points_1, desc_1] = sift(im1, 'Threshold', 0.015); % Do not change this threshold!
[points_2, desc_2] = sift(im2, 'Threshold', 0.015);
[points_3, desc_3] = sift(im3, 'Threshold', 0.015);
[points_4, desc_4] = sift(im4, 'Threshold', 0.015);

%% ToDo:

% Take image im1 as reference image (image 1) and compute the fundamental 
% matrices needed for computing the trajectory of point idx_car_I1
% (use the SIFT keypoints previously computed)

th = 2.0;

% im2 wrt im1
matches_im2_im1 = siftmatch(desc_2, desc_1);
p1 = [points_1(1:2, matches_im2_im1(2,:)); ones(1, length(matches_im2_im1))];
p2 = [points_2(1:2, matches_im2_im1(1,:)); ones(1, length(matches_im2_im1))];
[F_im2_im1, ~] = ransac_fundamental_matrix(p2, p1, th); 

% im3 wrt im1
matches_im3_im1 = siftmatch(desc_3, desc_1);
p1 = [points_1(1:2, matches_im3_im1(2,:)); ones(1, length(matches_im3_im1))];
p3 = [points_3(1:2, matches_im3_im1(1,:)); ones(1, length(matches_im3_im1))];
[F_im3_im1, ~] = ransac_fundamental_matrix(p3, p1, th); 

% im4 wrt im1
matches_im4_im1 = siftmatch(desc_4, desc_1);
p1 = [points_1(1:2, matches_im4_im1(2,:)); ones(1, length(matches_im4_im1))];
p4 = [points_4(1:2, matches_im4_im1(1,:)); ones(1, length(matches_im4_im1))];
[F_im4_im1, ~] = ransac_fundamental_matrix(p4, p1, th); 


%% Plot the car trajectory (keypoint idx_car_I1 in image 1)

% ToDo: complete the code

idx_car_I1 = 1197;
[~, index] = ismember(idx_car_I1, matches_im2_im1(2,:))
idx_car_I2 = matches_im2_im1(1, index); % ToDo: identify the corresponding point of idx_car_I1 in image 2
[~, index] = ismember(idx_car_I1, matches_im3_im1(2,:))
idx_car_I3 = matches_im3_im1(1, index); % ToDo: identify the corresponding point of idx_car_I1 in image 3
[~, index] = ismember(idx_car_I1, matches_im4_im1(2,:))
idx_car_I4 = matches_im4_im1(1, index); % ToDo: identify the corresponding point of idx_car_I1 in image 4

% coordinates (in image 1) of the keypoint idx_car_I1 (point in a van). 
% point1_1 is the projection of a 3D point in the 3D trajectory of the van
point1_1 = compute_point(points_1, idx_car_I1);
% coordinates (in image 1) of another 3D point in the same 3D trajectory of
% the van
point1_2 = [334 697 1]'; % (this is a given data)

% l1 is the projection of the 3D trajectory of keypoint idx_car_I1
% (it is the line that joins point1_1 and point1_2)
l1 = cross(point1_1, point1_2); % ToDo: compute the line
% plot the line
figure;imshow(im1);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(points_1(1,1197), points_1(2,1197), 'y*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 2
point2 = compute_point(points_2, idx_car_I2); %

% % ToDo: compute the epipolar line of point2 in the reference image
l2 = F_im2_im1*point2; %
% plot the epipolar line
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'c');
% ToDo: compute the projection of point idx_car_I2 in the reference image 
pi2 = cross(l1, l2); %
% plot this point
plot(pi2(1)/pi2(3), pi2(2)/pi2(3), 'c*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 3
point3 = compute_point(points_3, idx_car_I3); %
% ToDo: compute the epipolar line of point3 in the reference image
l3 = F_im3_im1*point3; %
% plot the epipolar line
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'b');
% ToDo: compute the projection of point idx_car_I3 in the reference image
pi3 = cross(l1, l3); %
plot(pi3(1)/pi3(3), pi3(2)/pi3(3), 'b*');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I1 in image 4
point4 = compute_point(points_4, idx_car_I4); %
% ToDo: compute the epipolar line of point4 in the reference image
l4 = F_im4_im1*point4; %
% plot the epipolar line
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'g');
% ToDo: compute the projection of point idx_car_I4 in the reference image
pi4 = cross(l1, l4); %
plot(pi4(1)/pi4(3), pi4(2)/pi4(3), 'g*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. OPTIONAL: Photo-sequencing with your own images

% 4.1 Take a set of images of a moving scene from different viewpoints at 
%     different time instants. At least two images have to be taken from
%     roughly the same location by the same camera.
%
% 4.2 Implement the first part (until line 16) of the Algorithm 1 of the 
%     Photo-sequencing paper with a selection of the detected dynamic
%     features. You may reuse the code generated for the previous question.
%
