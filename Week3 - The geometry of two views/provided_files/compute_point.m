function [ op ] = compute_point( p, idx )
% op is the projection of a 3D point in the 3D trajectory of the van
    op = [p(1:2, idx)' 1]';
end

