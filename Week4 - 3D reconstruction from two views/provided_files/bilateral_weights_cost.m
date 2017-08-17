function cost = bilateral_weights_cost(I1, I2)

center = ceil(size(I1)/2);

if(mod(size(I1,1),2) == 0)
    [p, q] = meshgrid(-center+1:center,-center+1:center);
    distance = sqrt(p.^2+q.^2);
else
    [p, q] = meshgrid(-center+1:center-1,-center+1:center-1);
    distance = sqrt(p.^2+q.^2);
end

gammac = 5;
gammap = 17.5;
T = 13.33;


num = sum(sum(exp(-abs(I1-I1(center(1),center(2)))/gammac-distance/gammap)...
    .*exp(-abs(I2-I2(center(1),center(2)))/gammac-distance/gammap)...
    .*min(abs(I1-I2),T)));
den = sum(sum(exp(-abs(I1-I1(center(1),center(2)))/gammac-distance/gammap)...
    .*exp(-abs(I2-I2(center(1),center(2)))/gammac-distance/gammap)));

cost = num/den;
  
end

% function [ cost ] = bilateral_weights_cost( w1, w2 )
% 
%     % gammac related to the strength of grouping by color similarity
%     % gammap is determined according to the size of the supportwindow 
%     % as gammap is proportional to window size
%     % T is the truncation value that controls the limit of the matching cost. 
%     gammac = 5;
%     gammap = 17.5;
%     T = 13.33;
%     
%     centroid = ceil(size(w1)/2);
% 
%     if(mod(size(w1,1),2) == 0)
%         [p, q] = meshgrid(-centroid+1:centroid,-centroid+1:centroid);
%         d = sqrt(p.^2+q.^2);
%     else
%         [p, q] = meshgrid(-centroid+1:centroid-1,-centroid+1:centroid-1);
%         d = sqrt(p.^2+q.^2);
%     end
%     
%     % Compute Support-Weight Based on the Strength of Grouping
%     Wpq_w1 = exp(-abs(w1-w1(centroid(1),centroid(2)))/gammac-d/gammap);
%     Wpq_w2 = exp(-abs(w2-w2(centroid(1),centroid(2)))/gammac-d/gammap);
%     e = min(abs(w1-w2),T);
%     cost = (sum(sum(Wpq_w1.*Wpq_w2.*e)) / sum(sum(Wpq_w1.*Wpq_w2)) );  
% end


