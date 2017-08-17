function [ outI ] = apply_H(I, H, reference)
% apply_H function that applies transformation H to image I
% returns image outI. 

sizeI = size(I);
sizeH = size(H);
[nrows, ncols, nchan] = size(I);

% Check dimensions 
if ((sizeH(1) ~= 3) || (sizeH(2) ~= 3))
    display('Matrix H has incorrect dimensions')
    return
end

% Fill in unset optional values.
switch nargin
    case 2
        reference = 'displayAllImage';
end

if strcmp(reference, 'displayAllImage')
    % Find corners in homogeneous coord.
    c1 = [1 1 1]';
    c2 = [ncols 1 1]';
    c3 = [1 nrows 1]';
    c4 = [ncols nrows 1]';

    Hc1=H*c1;
    Hc2=H*c2;
    Hc3=H*c3;
    Hc4=H*c4;

    Hc1 = Hc1/Hc1(3);
    Hc2 = Hc2/Hc2(3);
    Hc3 = Hc3/Hc3(3);
    Hc4 = Hc4/Hc4(3);

    xmin = round(min([Hc1(1) Hc2(1) Hc3(1) Hc4(1)]));
    xmax = round(max([Hc1(1) Hc2(1) Hc3(1) Hc4(1)]));
    ymin = round(min([Hc1(2) Hc2(2) Hc3(2) Hc4(2)]));
    ymax = round(max([Hc1(2) Hc2(2) Hc3(2) Hc4(2)]));
    
elseif strcmp(reference, 'keepOriginalPositions')
    xmin = 1;
    xmax = ncols;
    ymin = 1;
    ymax = nrows;
end

[X,Y] = meshgrid(xmin:xmax, ymin:ymax);
Hncols = xmax - xmin + 1;
Hnrows = ymax - ymin + 1;
Z = ones(Hnrows,Hncols);

XYZs = [X(:) Y(:) Z(:)]';

Hi = inv(H); 
HiXYZs = Hi * XYZs;
HX = reshape(HiXYZs(1,:), Hnrows, Hncols);
HY = reshape(HiXYZs(2,:), Hnrows, Hncols);
HZ = reshape(HiXYZs(3,:), Hnrows, Hncols);
HX = HX ./ HZ;
HY = HY ./ HZ;

outI = zeros(Hnrows,Hncols, nchan);
for c=1:nchan,
    outI(:,:,c) = interp2(double(I(:,:,c)), HX, HY, 'linear', 0);
end

% figure; imshow(uint8(I)); title('original image');
% figure; imshow(uint8(outI)); title('transformed image');

end
