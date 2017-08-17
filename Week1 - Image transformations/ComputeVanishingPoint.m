function [ point ] = ComputeVanishingPoint( line1, line2 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    point = cross(line1, line2);
    
    point = point / point(3);
end

