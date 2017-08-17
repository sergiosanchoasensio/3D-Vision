%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 1: Image rectification


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Applying image transformations

% ToDo: create the function  "apply_H" that gets as input a homography and
% an image and returns the image transformed by the homography.
% At some point you will need to interpolate the image values at some points,
% you may use the Matlab function "interp2" for that.

clear all;
close all;
clc;
%% 1.1. Similarities
I=imread('Data/0000_s.png'); % we have to be in the proper folder
figure(1); imshow(I); title('Original Image');

% ToDo: generate a matrix H which produces a similarity transformation
% Matrix H with translation
H = [1 0 100; 
     0 1 150;
     0 0  1 ];

I2 = apply_H(I, H, 'keepOriginalPositions');
figure(2); imshow(uint8(I2)); title('Planar transformations: Translation');

% Matrix H with rotation of 20º
H = [cosd(20) -sind(20) 0; 
     sind(20) cosd(20) 0;
     0 0 1];

I2 = apply_H(I, H);
figure(3); imshow(uint8(I2)); title('Planar transformations: Rotation');

% Matrix H with translation & rotation
H = [cosd(20) -sind(20) 100; 
     sind(20) cosd(20) 150;
     0 0 1];

I2 = apply_H(I, H, 'keepOriginalPositions');
figure(4); imshow(uint8(I2)); title('Planar transformations: Translation & Rotation');

% Matrix H with scale
H = [2 0 0; 
     0 2 0;
     0 0 1];

I2 = apply_H(I, H);
figure(5); imshow(uint8(I2)); title('Planar transformations: Scale');


%% 1.2. Affinities

% ToDo: generate a matrix H which produces an affine transformation
H = [1.15  -0.3  160; 
     0.3   0.65  -65;
     0      0     1];

I2 = apply_H(I, H);
figure(6); imshow(uint8(I2)); title('Planar transformations: Affine')

% ToDo: decompose the affinity in four transformations: two
% rotations, a scale, and a translation
A = H(1:size(H,1)-1,1:size(H,2)-1);
[U,S,V] = svd(A); %performs a singular value decomposition of matrix A, such that A = U*S*V'.
R_theta = U*V';
R_phi = V';
D = S;
t = [H(1,size(H,1)); H(2,size(H,1))];

% ToDo: verify that the product of the four previous transformations
% produces the same matrix H as above

% Compute A, following Lecture 1, Page 43 of 48
A = R_theta*R_phi'*D*R_phi;
% Add the translation terms to H
H_new = [A t; 0 0 1];
I2 = apply_H(I, H_new);
figure(7); imshow(uint8(I2)); title('Planar transformations: Affine (From SVD decomposition)')

% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before
Identity = [1 0; 0 1];
H_R_phi = [R_phi [0; 0]; 0 0 1];
H_D = [D [0; 0]; 0 0 1];
H_R_theta = [R_theta [0; 0]; 0 0 1];
H_t = [Identity t; 0 0 1];

I_aux = apply_H(I, H_R_theta);
I_aux = apply_H(I_aux, H_R_phi');
I_aux = apply_H(I_aux, H_D);
I_aux = apply_H(I_aux, H_R_phi);
I_aux = apply_H(I_aux, H_t);
figure(8); imshow(uint8(I_aux)); title('Planar transformations: Affine (Apply matrices independently)')

%% 1.3 Projective transformations (homographies)
% ToDo: generate a matrix H which produces a projective transformation
H = [0.7 -1 -0.5; 0.2 0.8 -1; 0.00004 0 -1];
I2 = apply_H(I, H);
figure(9); imshow(uint8(I2));title('Projective transformations')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification

% choose the image points
I = imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
i = 424;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 240;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 712;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 565;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';
i = 227;
p9 = [A(i,1) A(i,2) 1]';
p10 = [A(i,3) A(i,4) 1]';
i = 576;
p11 = [A(i,1) A(i,2) 1]';
p12 = [A(i,3) A(i,4) 1]';
i = 534;
p13 = [A(i,1) A(i,2) 1]';
p14 = [A(i,3) A(i,4) 1]';

% ToDo: compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = computeLine( p1, p2);
l2 = computeLine( p3, p4);
l3 = computeLine( p5, p6);
l4 = computeLine( p7, p8);

% show the chosen lines in the image
figure(10);imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');
% plot(t, -(l5(1)*t + l5(3)) / l5(2), 'y');
% plot(t, -(l6(1)*t + l6(3)) / l6(2), 'y');

% ToDo: compute the homography that affinely rectifies the image
% Compute vanishing points
vp_12 = cross(l1, l2);
vp_12 = vp_12 / vp_12(3);
vp_34 = cross(l3, l4);
vp_34 = vp_34 / vp_34(3);
l_inf = cross(vp_12,vp_34);
l_inf = l_inf / l_inf(3);

H_ap = [1 0 0; 0 1 0; l_inf];
I_ap = apply_H(I, H_ap);
figure(11); imshow(uint8(I_ap)); title('Affine rectification via the vanishing line')

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
points = [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14];
newPoints  = apply_H_toPoints( H_ap, points );

lr1 = computeLine( newPoints(:,1),  newPoints(:,2));
lr2 = computeLine( newPoints(:,3),  newPoints(:,4));
lr3 = computeLine( newPoints(:,5),  newPoints(:,6));
lr4 = computeLine( newPoints(:,7),  newPoints(:,8));
lr5 = computeLine( newPoints(:,9),  newPoints(:,12));
lr6 = computeLine( newPoints(:,10), newPoints(:,13));

% show the transformed lines in the transformed image
figure(12);imshow(uint8(I_ap)); title('Affine rectification. Rectificated lines')
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation

% get normalized version of lines (two coordinates)
norml1 = [l1(1)/l1(3), l1(2)/l1(3)];
norml2 = [l2(1)/l2(3), l2(2)/l2(3)];
norml3 = [l3(1)/l1(3), l3(2)/l3(3)];
norml4 = [l4(1)/l4(3), l4(2)/l4(3)];

normlr1 = [lr1(1)/lr1(3), lr1(2)/lr1(3)];
normlr2 = [lr2(1)/lr2(3), lr2(2)/lr2(3)];
normlr3 = [lr3(1)/lr1(3), lr3(2)/lr3(3)];
normlr4 = [lr4(1)/lr4(3), lr4(2)/lr4(3)];

angle13 = acosd(dot(norml1,norml3)/(norm(norml1)*norm(norml3)));
angler13 = acosd(dot(normlr1,normlr3)/(norm(normlr1)*norm(normlr3)));

angle14 = acosd(dot(norml1,norml4)/(norm(norml1)*norm(norml4)));
angler14 = acosd(dot(normlr1,normlr4)/(norm(normlr1)*norm(normlr4)));

angle23= acosd(dot(norml2,norml3)/(norm(norml2)*norm(norml3)));
angler23 = acosd(dot(normlr2,normlr3)/(norm(normlr2)*norm(normlr3)));

angle24= acosd(dot(norml2,norml4)/(norm(norml2)*norm(norml4)));
angler24 = acosd(dot(normlr2,normlr4)/(norm(normlr2)*norm(normlr4)));

disp(['Upper left corner before transformation: ' , num2str(angle13), 'º']);
disp(['Upper left corner after transformation: ' , num2str(angler13), 'º']);
disp(' ');

disp(['Upper right corner before transformation: ' , num2str(angle14), 'º']);
disp(['Upper right corner after transformation: ' , num2str(angler14), 'º']);
disp(' ');

disp(['Lower left corner before transformation: ' , num2str(angle23), 'º']);
disp(['Lower left corner after transformation: ' , num2str(angler23), 'º']);
disp(' ');

disp(['Lower right corner before transformation: ' , num2str(angle24), 'º']);
disp(['Lower right corner after transformation: ' , num2str(angler24), 'º']);
disp(' ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.
lr1 = computeLine( newPoints(:,1),  newPoints(:,2));
lr3 = computeLine( newPoints(:,5),  newPoints(:,6));
lr5 = computeLine( newPoints(:,9) , newPoints(:,12));
lr6 = computeLine( newPoints(:,10), newPoints(:,13));

% normalize orthogonal lines so that they have 2 elements
lr1 = lr1 / lr1(3);
lr3 = lr3 / lr3(3);
lr5 = lr5 / lr5(3);
lr6 = lr6 / lr6(3);

figure(13);imshow(uint8(I_ap)); title('Before metric rectification.')
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr5(1)*t + lr5(3)) / lr5(2), 'r');
plot(t, -(lr6(1)*t + lr6(3)) / lr6(2), 'r');

n = [lr5(1)*lr6(1), lr5(1)*lr6(2) + lr5(2)*lr6(1), lr5(2)*lr6(2) ; ...
     lr1(1)*lr3(1), lr1(1)*lr3(2) + lr1(2)*lr3(1), lr1(2)*lr3(2)];
 
s = null(n);
S = [s(1) s(2);s(2) s(3)];
% Cholesky decomposition
K = chol(S,'upper');

H_sa = [K [0; 0]; 0 0 1];
H_sa = inv(H_sa);
I_sa = apply_H(I_ap, H_sa);
figure(14); imshow(uint8(I_sa)); title('Metric rectification via orthogonal lines.')

% Compute image with lines
newPoints_sa  = apply_H_toPoints( H_sa, newPoints );

lrr1 = computeLine( newPoints_sa(:,1),  newPoints_sa(:,2));
lrr3 = computeLine( newPoints_sa(:,5),  newPoints_sa(:,6));
lrr5 = computeLine( newPoints_sa(:,9) , newPoints_sa(:,12));
lrr6 = computeLine( newPoints_sa(:,10), newPoints_sa(:,13));

% normalize orthogonal lines so that they have 2 elements
lrr1 = lrr1 / lrr1(3);
lrr3 = lrr3 / lrr3(3);
lrr5 = lrr5 / lrr5(3);
lrr6 = lrr6 / lrr6(3);

% show the transformed lines in the transformed image
figure(15);imshow(uint8(I_sa)); title('Metric rectification. Rectificated lines')
hold on;
t=1:0.1:1000;
plot(t, -(lrr1(1)*t + lrr1(3)) / lrr1(2), 'y');
plot(t, -(lrr3(1)*t + lrr3(3)) / lrr3(2), 'y');
plot(t, -(lrr5(1)*t + lrr5(3)) / lrr5(2), 'r');
plot(t, -(lrr6(1)*t + lrr6(3)) / lrr6(2), 'r');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation
% get normalized version of lines (two coordinates)
normlr1 = [lr1(1), lr1(2)];
normlr3 = [lr3(1), lr3(2)];
normlr5 = [lr5(1), lr5(2)];
normlr6 = [lr6(1), lr6(2)];

normlrr1 = [lrr1(1), lrr1(2)];
normlrr3 = [lrr3(1), lrr3(2)];
normlrr5 = [lrr5(1), lrr5(2)];
normlrr6 = [lrr6(1), lrr6(2)];

angle13 = acosd(dot(normlr1,normlr3)/(norm(normlr1)*norm(normlr3)));
angler13 = acosd(dot(normlrr1,normlrr3)/(norm(normlrr1)*norm(normlrr3)));

angle56 = acosd(dot(normlr5,normlr6)/(norm(normlr5)*norm(normlr6)));
angler56 = acosd(dot(normlrr5,normlrr6)/(norm(normlrr5)*norm(normlrr6)));

disp(['Crossing point yellow lines before transformation: ' , num2str(angle13), 'º']);
disp(['Crossing point yellow lines after transformation: ' , num2str(angler13), 'º']);
disp(' ');

disp(['Crossing point red lines before transformation: ' , num2str(angle56), 'º']);
disp(['Crossing point red lines after transformation: ' , num2str(angler56), 'º']);
disp(' ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)

% choose the image points
I = imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
i = 424;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 240;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 712;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 565;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';
i = 227;
p9 = [A(i,1) A(i,2) 1]';
p10 = [A(i,3) A(i,4) 1]';
i = 576;
p11 = [A(i,1) A(i,2) 1]';
p12 = [A(i,3) A(i,4) 1]';
i = 534;
p13 = [A(i,1) A(i,2) 1]';
p14 = [A(i,3) A(i,4) 1]';

% Construct 5 pairs of orthogonal lines
l1  = computeLine( p1  , p2);
l2  = computeLine( p5  , p6);

l3  = computeLine( p3  , p4);
l4  = computeLine( p7  , p8);

l5  = computeLine( p9  , p12);
l6  = computeLine( p10 , p13);

l7  = computeLine( p9  , p10);
l8  = computeLine( p13 , p14);

l9  = computeLine( p2  , p5);
l10 = computeLine( p1  , p8);

% normalize lines
l1  = l1 / l1(3);
l2  = l2 / l2(3);
l3  = l3 / l3(3);
l4  = l4 / l4(3);
l5  = l5 / l5(3);
l6  = l6 / l6(3);
l7  = l7 / l7(3);
l8  = l8 / l8(3);
l9  = l9 / l9(3);
l10 = l10 / l10(3);

% show the chosen lines in the image
figure(16);imshow(I);
hold on;
t=1:0.1:1000; title('Pairs of orthogonal lines.')
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'r');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'r');
plot(t, -(l5(1)*t + l5(3)) / l5(2), 'm');
plot(t, -(l6(1)*t + l6(3)) / l6(2), 'm');
plot(t, -(l7(1)*t + l7(3)) / l7(2), 'g');
plot(t, -(l8(1)*t + l8(3)) / l8(2), 'g');
plot(t, -(l9(1)*t + l9(3)) / l9(2), 'b');
plot(t, -(l10(1)*t + l10(3)) / l10(2), 'b');

% Find conic coefficient matrix C*
c = [l1(1)*l2(1) , (l1(1)*l2(2)  + l1(2)*l2(1))/2 , l1(2)*l2(2) , (l1(1)+l2(1))/2 , (l1(2)+l2(2))/2 , 1; ...
     l3(1)*l4(1) , (l3(1)*l4(2)  + l3(2)*l4(1))/2 , l3(2)*l4(2) , (l3(1)+l4(1))/2 , (l3(2)+l4(2))/2 , 1; ...
     l5(1)*l6(1) , (l5(1)*l6(2)  + l5(2)*l6(1))/2 , l5(2)*l6(2) , (l5(1)+l6(1))/2 , (l5(2)+l6(2))/2 , 1; ...
     l7(1)*l8(1) , (l7(1)*l8(2)  + l7(2)*l8(1))/2 , l7(2)*l8(2) , (l7(1)+l8(1))/2 , (l7(2)+l8(2))/2 , 1; ...
     l9(1)*l10(1), (l9(1)*l10(2) + l9(2)*l10(1))/2, l9(2)*l10(2), (l9(1)+l10(1))/2, (l9(2)+l10(2))/2, 1];
     
c = null(c);
% Get the "projected" conic matrix (C*') from c coefficients. (Cprojected)
Cprojected = [c(1) c(2)/2 c(4)/2 ; c(2)/2 c(3) c(5)/2; c(4)/2 c(5)/2 c(6)];

[U,S,V] = svd(Cprojected);
H = inv(U);

I_out = apply_H(I, H);
figure(17); imshow(uint8(I_out)); title('Direct metric rectification via orthogonal lines.')
%% 5. OPTIONAL: Affine Rectification of the left facade of image 0000
% choose the image points
I = imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
i = 493;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 186;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 48;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 508;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

i = 192; %6; % 62
p9 = [A(i,1) A(i,2) 1]';
p10 = [A(i,3) A(i,4) 1]';
i = 195; %7; % 67 
p11 = [A(i,1) A(i,2) 1]';
p12 = [A(i,3) A(i,4) 1]';

% compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = computeLine( p1, p2);
l2 = computeLine( p3, p4);
l3 = computeLine( p5, p6);
l4 = computeLine( p7, p8);
l5 = computeLine( p9, p10);
l6 = computeLine( p11, p12);

% show the chosen lines in the image
figure(17);imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');

plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% plot(t, -(l5(1)*t + l5(3)) / l5(2), 'r');
% plot(t, -(l6(1)*t + l6(3)) / l6(2), 'r');


% ToDo: compute the homography that affinely rectifies the image
% Compute vanishing points
vp_12 = cross(l1, l2);
vp_12 = vp_12 / vp_12(3);
vp_34 = cross(l3, l4);
vp_34 = vp_34 / vp_34(3);
l_inf = cross(vp_12,vp_34);
l_inf = l_inf / l_inf(3);

H_ap = [1 0 0; 0 1 0; l_inf];
I_ap = apply_H(I, H_ap, 'keepOriginalPositions');
figure(18); imshow(uint8(I_ap)); title('Affine rectification via the vanishing line')

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
points = [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12];
newPoints  = apply_H_toPoints( H_ap, points );

lr1 = computeLine( newPoints(:,1),  newPoints(:,2));
lr2 = computeLine( newPoints(:,3),  newPoints(:,4));
lr3 = computeLine( newPoints(:,5),  newPoints(:,6));
lr4 = computeLine( newPoints(:,7),  newPoints(:,8));
lr5 = computeLine( newPoints(:,9),  newPoints(:,10));
lr6 = computeLine( newPoints(:,11),  newPoints(:,12));

% show the transformed lines in the transformed image
figure(19);imshow(uint8(I_ap)); title('Affine rectification. Rectificated lines')
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');


% to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation

norml1 = [l1(1)/l1(3), l1(2)/l1(3)];
norml2 = [l2(1)/l2(3), l2(2)/l2(3)];
norml3 = [l3(1)/l1(3), l3(2)/l3(3)];
norml4 = [l4(1)/l4(3), l4(2)/l4(3)];

normlr1 = [lr1(1)/lr1(3), lr1(2)/lr1(3)];
normlr2 = [lr2(1)/lr2(3), lr2(2)/lr2(3)];
normlr3 = [lr3(1)/lr1(3), lr3(2)/lr3(3)];
normlr4 = [lr4(1)/lr4(3), lr4(2)/lr4(3)];

angle13 = acosd(dot(norml1,norml3)/(norm(norml1)*norm(norml3)));
angler13 = acosd(dot(normlr1,normlr3)/(norm(normlr1)*norm(normlr3)));

angle14 = acosd(dot(norml1,norml4)/(norm(norml1)*norm(norml4)));
angler14 = acosd(dot(normlr1,normlr4)/(norm(normlr1)*norm(normlr4)));

angle23 = acosd(dot(norml2,norml3)/(norm(norml2)*norm(norml3)));
angler23 = acosd(dot(normlr2,normlr3)/(norm(normlr2)*norm(normlr3)));

angle24 = acosd(dot(norml2,norml4)/(norm(norml2)*norm(norml4)));
angler24 = acosd(dot(normlr2,normlr4)/(norm(normlr2)*norm(normlr4)));

disp(['Upper left corner before transformation: ' , num2str(angle13), 'º']);
disp(['Upper left corner after transformation: ' , num2str(angler13), 'º']);
disp(' ');

disp(['Upper right corner before transformation: ' , num2str(angle14), 'º']);
disp(['Upper right corner after transformation: ' , num2str(angler14), 'º']);
disp(' ');

disp(['Lower left corner before transformation: ' , num2str(angle23), 'º']);
disp(['Lower left corner after transformation: ' , num2str(angler23), 'º']);
disp(' ');

disp(['Lower right corner before transformation: ' , num2str(angle24), 'º']);
disp(['Lower right corner after transformation: ' , num2str(angler24), 'º']);
disp(' ');


%% 6. OPTIONAL: Metric Rectification of the left facade of image 0000


lr1 = computeLine( newPoints(:,1),  newPoints(:,2));
lr3 = computeLine( newPoints(:,5),  newPoints(:,6));
lr5 = computeLine( newPoints(:,3), newPoints(:,2));
lr6 = computeLine( newPoints(:,4), newPoints(:,1));


% l1 = computeLine( p1, p2);
% l3 = computeLine( p5, p6);
% l5 = computeLine( p3, p2);
% l6 = computeLine( p4, p1);


figure(20);imshow(uint8(I_ap));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr5(1)*t + lr5(3)) / lr5(2), 'r');
plot(t, -(lr6(1)*t + lr6(3)) / lr6(2), 'r');

% normalize orthogonal lines so that they have 2 elements
lr1 = lr1 / lr1(3);
lr3 = lr3 / lr3(3);
lr5 = lr5 / lr5(3);
lr6 = lr6 / lr6(3);

figure(21);imshow(uint8(I_ap)); title('Before metric rectification.')
hold on;
t=1:0.1:1000;

plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr5(1)*t + lr5(3)) / lr5(2), 'r');
plot(t, -(lr6(1)*t + lr6(3)) / lr6(2), 'r');

n = [lr5(1)*lr6(1), lr5(1)*lr6(2) + lr5(2)*lr6(1), lr5(2)*lr6(2) ; ...
     lr1(1)*lr3(1), lr1(1)*lr3(2) + lr1(2)*lr3(1), lr1(2)*lr3(2)];
 
s = null(n);
S = [s(1) s(2);s(2) s(3)];
% Cholesky decomposition
K = chol(S,'upper');


H_sa = [K [0; 0]; 0 0 1];
H_sa = inv(H_sa);
I_sa = apply_H(I_ap, H_sa);
figure(21); imshow(uint8(I_sa)); title('Metric rectification via orthogonal lines.')

% Compute image with lines
newPoints_sa  = apply_H_toPoints( H_sa, newPoints );

lrr1 = computeLine( newPoints_sa(:,1),  newPoints_sa(:,2));
lrr3 = computeLine( newPoints_sa(:,5),  newPoints_sa(:,6));
lrr5 = computeLine( newPoints_sa(:,3) , newPoints_sa(:,2));
lrr6 = computeLine( newPoints_sa(:,4), newPoints_sa(:,1));

% normalize orthogonal lines so that they have 2 elements
lrr1 = lrr1 / lrr1(3);
lrr3 = lrr3 / lrr3(3);
lrr5 = lrr5 / lrr5(3);
lrr6 = lrr6 / lrr6(3);

% show the transformed lines in the transformed image
figure(22);imshow(uint8(I_sa)); title('Metric rectification. Rectificated lines')
hold on;
t=1:0.1:1000;
plot(t, -(lrr1(1)*t + lrr1(3)) / lrr1(2), 'y');
plot(t, -(lrr3(1)*t + lrr3(3)) / lrr3(2), 'y');
plot(t, -(lrr5(1)*t + lrr5(3)) / lrr5(2), 'r');
plot(t, -(lrr6(1)*t + lrr6(3)) / lrr6(2), 'r');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation
% get normalized version of lines (two coordinates)
normlr1 = [lr1(1), lr1(2)];
normlr3 = [lr3(1), lr3(2)];
normlr5 = [lr5(1), lr5(2)];
normlr6 = [lr6(1), lr6(2)];

normlrr1 = [lrr1(1), lrr1(2)];
normlrr3 = [lrr3(1), lrr3(2)];
normlrr5 = [lrr5(1), lrr5(2)];
normlrr6 = [lrr6(1), lrr6(2)];

angle13 = acosd(dot(normlr1,normlr3)/(norm(normlr1)*norm(normlr3)));
angler13 = acosd(dot(normlrr1,normlrr3)/(norm(normlrr1)*norm(normlrr3)));

angle56 = acosd(dot(normlr5,normlr6)/(norm(normlr5)*norm(normlr6)));
angler56 = acosd(dot(normlrr5,normlrr6)/(norm(normlrr5)*norm(normlrr6)));

disp(['Crossing point yellow lines before transformation: ' , num2str(angle13), 'º']);
disp(['Crossing point yellow lines after transformation: ' , num2str(angler13), 'º']);
disp(' ');

disp(['Crossing point red lines before transformation: ' , num2str(angle56), 'º']);
disp(['Crossing point red lines after transformation: ' , num2str(angler56), 'º']);
disp(' ');

%% 7. OPTIONAL: Affine Rectification of the left facade of image 0001
I = imread('Data/0001_s.png');
A = load('Data/0001_s_info_lines.txt');

% indices of lines
i = 614;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 159;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 645;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 541;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

i = 1;
p9 = [A(i,1) A(i,2) 1]';
p10 = [A(i,3) A(i,4) 1]';
i = 5;
p11 = [A(i,1) A(i,2) 1]';
p12 = [A(i,3) A(i,4) 1]';

% compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = computeLine( p1, p2);
l2 = computeLine( p3, p4);
l3 = computeLine( p5, p6);
l4 = computeLine( p7, p8);
l5 = computeLine( p9, p10);
l6 = computeLine( p11, p12);

% show the chosen lines in the image
figure(23);imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'y');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'y');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'y');

% plot(t, -(l5(1)*t + l5(3)) / l5(2), 'r');
% plot(t, -(l6(1)*t + l6(3)) / l6(2), 'r');

% ToDo: compute the homography that affinely rectifies the image
% Compute vanishing points
vp_12 = cross(l1, l2);
vp_12 = vp_12 / vp_12(3);
vp_34 = cross(l3, l4);
vp_34 = vp_34 / vp_34(3);
l_inf = cross(vp_12,vp_34);
l_inf = l_inf / l_inf(3);

H_ap = [1 0 0; 0 1 0; l_inf];
I_ap = apply_H(I, H_ap, 'keepOriginalPositions');
figure(24); imshow(uint8(I_ap)); title('Affine rectification via the vanishing line')

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4
points = [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12];
newPoints  = apply_H_toPoints( H_ap, points );

lr1 = computeLine( newPoints(:,1),  newPoints(:,2));
lr2 = computeLine( newPoints(:,3),  newPoints(:,4));
lr3 = computeLine( newPoints(:,5),  newPoints(:,6));
lr4 = computeLine( newPoints(:,7),  newPoints(:,8));
lr5 = computeLine( newPoints(:,9),  newPoints(:,10));
lr6 = computeLine( newPoints(:,11),  newPoints(:,12));

% show the transformed lines in the transformed image
figure(25);imshow(uint8(I_ap)); title('Affine rectification. Rectificated lines')
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'y');

% to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation

norml1 = [l1(1)/l1(3), l1(2)/l1(3)];
norml2 = [l2(1)/l2(3), l2(2)/l2(3)];
norml3 = [l3(1)/l1(3), l3(2)/l3(3)];
norml4 = [l4(1)/l4(3), l4(2)/l4(3)];

normlr1 = [lr1(1)/lr1(3), lr1(2)/lr1(3)];
normlr2 = [lr2(1)/lr2(3), lr2(2)/lr2(3)];
normlr3 = [lr3(1)/lr1(3), lr3(2)/lr3(3)];
normlr4 = [lr4(1)/lr4(3), lr4(2)/lr4(3)];

angle13 = acosd(dot(norml1,norml3)/(norm(norml1)*norm(norml3)));
angler13 = acosd(dot(normlr1,normlr3)/(norm(normlr1)*norm(normlr3)));

angle14 = acosd(dot(norml1,norml4)/(norm(norml1)*norm(norml4)));
angler14 = acosd(dot(normlr1,normlr4)/(norm(normlr1)*norm(normlr4)));

angle23 = acosd(dot(norml2,norml3)/(norm(norml2)*norm(norml3)));
angler23 = acosd(dot(normlr2,normlr3)/(norm(normlr2)*norm(normlr3)));

angle24 = acosd(dot(norml2,norml4)/(norm(norml2)*norm(norml4)));
angler24 = acosd(dot(normlr2,normlr4)/(norm(normlr2)*norm(normlr4)));

disp(['Upper left corner before transformation: ' , num2str(angle13), 'º']);
disp(['Upper left corner after transformation: ' , num2str(angler13), 'º']);
disp(' ');

disp(['Upper right corner before transformation: ' , num2str(angle14), 'º']);
disp(['Upper right corner after transformation: ' , num2str(angler14), 'º']);
disp(' ');

disp(['Lower left corner before transformation: ' , num2str(angle23), 'º']);
disp(['Lower left corner after transformation: ' , num2str(angler23), 'º']);
disp(' ');

disp(['Lower right corner before transformation: ' , num2str(angle24), 'º']);
disp(['Lower right corner after transformation: ' , num2str(angler24), 'º']);
disp(' ');

%% 8. OPTIONAL: Metric Rectification of the left facade of image 0001

lr1 = computeLine(newPoints(:,1), newPoints(:,2));
lr3 = computeLine(newPoints(:,3), newPoints(:,1));
lr5 = computeLine(newPoints(:,3), newPoints(:,2));
lr6 = computeLine(newPoints(:,4), newPoints(:,1));

figure(26);imshow(uint8(I_ap));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr5(1)*t + lr5(3)) / lr5(2), 'r');
plot(t, -(lr6(1)*t + lr6(3)) / lr6(2), 'r');

% normalize orthogonal lines so that they have 2 elements
lr1 = lr1 / lr1(3);
lr3 = lr3 / lr3(3);
lr5 = lr5 / lr5(3);
lr6 = lr6 / lr6(3);

figure(27);imshow(uint8(I_ap)); title('Before metric rectification.')
hold on;
t=1:0.1:1000;

plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'y');
plot(t, -(lr5(1)*t + lr5(3)) / lr5(2), 'r');
plot(t, -(lr6(1)*t + lr6(3)) / lr6(2), 'r');

n = [lr5(1)*lr6(1), lr5(1)*lr6(2) + lr5(2)*lr6(1), lr5(2)*lr6(2) ; ...
     lr1(1)*lr3(1), lr1(1)*lr3(2) + lr1(2)*lr3(1), lr1(2)*lr3(2)];
 
s = null(n);
S = [s(1) s(2);s(2) s(3)];
% Cholesky decomposition
K = chol(S,'upper');


H_sa = [K [0; 0]; 0 0 1];
H_sa = inv(H_sa);
I_sa = apply_H(I_ap, H_sa);
figure(28); imshow(uint8(I_sa)); title('Metric rectification via orthogonal lines.')

% Compute image with lines
newPoints_sa  = apply_H_toPoints( H_sa, newPoints );

lrr1 = computeLine( newPoints_sa(:,1),  newPoints_sa(:,2));
lrr3 = computeLine( newPoints_sa(:,1),  newPoints_sa(:,3));
lrr5 = computeLine( newPoints_sa(:,2) , newPoints_sa(:,3));
lrr6 = computeLine( newPoints_sa(:,1), newPoints_sa(:,4));

% normalize orthogonal lines so that they have 2 elements
lrr1 = lrr1 / lrr1(3);
lrr3 = lrr3 / lrr3(3);
lrr5 = lrr5 / lrr5(3);
lrr6 = lrr6 / lrr6(3);

% show the transformed lines in the transformed image
figure(29);imshow(uint8(I_sa)); title('Metric rectification. Rectificated lines')
hold on;
t=1:0.1:10000;
plot(t, -(lrr1(1)*t + lrr1(3)) / lrr1(2), 'y');
plot(t, -(lrr3(1)*t + lrr3(3)) / lrr3(2), 'y');
plot(t, -(lrr5(1)*t + lrr5(3)) / lrr5(2), 'r');
plot(t, -(lrr6(1)*t + lrr6(3)) / lrr6(2), 'r');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation
% get normalized version of lines (two coordinates)
normlr1 = [lr1(1), lr1(2)];
normlr3 = [lr3(1), lr3(2)];
normlr5 = [lr5(1), lr5(2)];
normlr6 = [lr6(1), lr6(2)];

normlrr1 = [lrr1(1), lrr1(2)];
normlrr3 = [lrr3(1), lrr3(2)];
normlrr5 = [lrr5(1), lrr5(2)];
normlrr6 = [lrr6(1), lrr6(2)];

angle13 = acosd(dot(normlr1,normlr3)/(norm(normlr1)*norm(normlr3)));
angler13 = acosd(dot(normlrr1,normlrr3)/(norm(normlrr1)*norm(normlrr3)));

angle56 = acosd(dot(normlr5,normlr6)/(norm(normlr5)*norm(normlr6)));
angler56 = acosd(dot(normlrr5,normlrr6)/(norm(normlrr5)*norm(normlrr6)));

disp(['Crossing point yellow lines before transformation: ' , num2str(angle13), 'º']);
disp(['Crossing point yellow lines after transformation: ' , num2str(angler13), 'º']);
disp(' ');

disp(['Crossing point red lines before transformation: ' , num2str(angle56), 'º']);
disp(['Crossing point red lines after transformation: ' , num2str(angler56), 'º']);
disp(' ');