function mvp_all = refine_detections3(mvp_all, lines_lsd, params);
% refines VP detections using lines from LSD
D = size(mvp_all,2);

mvp_refined = zeros(D,2);
for i =1:D
    vp = mvp_all(:,i)';
    vp = refine_vp(lines_lsd, vp, params);
    mvp_refined(i,:) = vp;
end

mvp_all = mvp_refined';
end

%% iterator function
function  vp = refine_vp( lines,  vp, params);

% given a cluster of line segments, and two segments indicated by p1 and
% p2, obtain the main vanishing point determined by the segments

THRESHOLD = params.REFINE_THRESHOLD;
% THRESHOLD = THRESHOLD-(THRESHOLD/1.5)

H = params.H;
W = params.W;

vp_orig = vp;
vp_tmp = vp;

[vp] = refine_vp_iteration(lines,  vp, THRESHOLD,H, W);


variation = norm(vp-vp_orig)/norm(vp_orig);
fprintf('variation: %f\n', variation);
if (variation > params.VARIATION_THRESHOLD)
    fprintf('vp changed too much (%f)... probably unstable conditions, keeping initial vp\n', variation)
    vp = vp_orig;
end

end

%% iteration function
function [vp] = refine_vp_iteration(lines,  vp, THRESHOLD, H,  W)
% finds intersection of each line in cluster with vanishing point segments

%debug
% if threshold finishes in .1 use all lines
z=1:length(lines);

vp_orig = vp;


z2= [];

X = [0 W];

mp =[lines(:,1)+(lines(:,3)-lines(:,1))/2, lines(:,2)+(lines(:,4)-lines(:,2))/2];

L = size(lines,1);
O = ones(L,1);
Z = zeros(L,1);
vpmat = my_repmat2(vp,[L 1]);


VP = my_cross([mp O], [vpmat O]);
VP3 = my_repmat(VP(:,3),[1 3]);
VP = VP./VP3;

mp_vp = [VP(:,1) VP(:,2)];


a = VP(:,1);
b = VP(:,2);

% get angle betwen lines
lt2 = [Z -1./b W*O -W*a./b-1./b];

A = lines(:,3:4)-lines(:,1:2);
B = lt2(:,3:4)-lt2(:,1:2);

normA = sqrt(A(:,1).^2+A(:,2).^2);
normB = sqrt(B(:,1).^2+B(:,2).^2);

A = A./my_repmat(normA,[1 2]);
B = B./my_repmat(normB,[1 2]);

angle = acos(dot(A',B')');
angle = real(angle); % numerical errors
angle = min(angle, pi-angle);
    
    angle = angle*180/pi;
    
    z2 = find(angle<THRESHOLD);

    Z = length(z);

% for j=1:Z
%     i=z(j);
%     lt = lines(i,:);
%     
%     % get line from segment middle point to vp
%     mp =[lt(1)+(lt(3)-lt(1))/2, lt(2)+(lt(4)-lt(2))/2];
%     
%     
%     mp_vp = get_intersection([mp 1], [vp 1]);
%     
%     a = mp_vp(1);
%     b = mp_vp(2);
%     
%     % get angle betwen lines
%     lt2 = [0 -1/b W -W*a/b-1/b];
%     
%     A = lt(3:4)-lt(1:2);
%     B = lt2(3:4)-lt2(1:2);
%     
%     
%     angle = acos(dot(A,B)/norm(A)/norm(B));
%     angle = min(angle, pi-angle);
%     
%     angle = angle*180/pi;
%     
%     if angle< THRESHOLD
%         z2 = [z2 i];
%     end
%     
%     
%     
% end

% fprintf('found %i aligned segments out of %i\n', length(z2), Z);

%% obtain a refined VP estimate from sub-cluster z2
lengths = sqrt(sum(([lines(:,1) lines(:,2)] - [lines(:,3) lines(:,4)]).^2,2));
weights = lengths/max(lengths);
lis=line_to_homogeneous(lines);

Z2 = length(z2);

Q=zeros(3);
Is = [1 0 0; 0 1 0; 0 0 0];

% for I=1:Z2
%     i = z2(I);
%     li = lis(i,:)'; %line_to_homogeneous(lines(i,:))';
%     q = weights(i)^2*(li*li')/(li'*Is*li);
%     Q = Q+q;
% end

% 
l2 = lis(z2,:)';
w2 = weights(z2).^2;
w2 = my_repmat(w2,[1 3])';

b = dot((l2'*Is)',l2)'; %diag(l2'*Is*l2);
b = my_repmat(b,[1 3])';
Q = (l2./b.*w2)*l2';


p = [0 0 1]';

A = [2*Q -p];

vp = null(A);

% fprintf('cond(A) = %2.2e, cond(A^2) = %2.2e, cond(Q) = %2.2e, vp(3) = %2.2e\n', cond(A), cond(A'*A), cond(Q),vp(3,1));

vp = vp(1:2,1)/vp(3,1);
vp = vp';
%
%
% if    cond(A)>1e9 %cond(A'*A)>1e22 % &
%     
%     fprintf('matrix ill conditioned, probably mostly parallel lines, obtain VP by averaging...\n');
%     vp = vp_orig;
% end


end



