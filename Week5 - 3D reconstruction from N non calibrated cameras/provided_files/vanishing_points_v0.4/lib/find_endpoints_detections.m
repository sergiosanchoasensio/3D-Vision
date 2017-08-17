function detections = find_endpoints_detections(lines,z,params)
% find detections of alignments among line segment endpoints

lines2 = lines(z,:);

H = params.H;
W = params.W;

line_size = params.line_size;

points_segments = [[lines2(:,1); lines2(:,3)], H-[lines2(:,2); lines2(:,4)]];

params.endpoints = 1;
[detections, m, b] = find_detections(points_segments, params);
params.endpoints = 0;