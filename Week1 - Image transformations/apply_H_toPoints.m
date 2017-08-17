function [ transformedPoints ] = apply_H_toPoints( H, pointsArray )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    transformedPoints = [];
    for i = 1 : size(pointsArray,2)
        transformedPoints  = [transformedPoints computeTransformedPoints(H, pointsArray(:,i))];
    end
end

