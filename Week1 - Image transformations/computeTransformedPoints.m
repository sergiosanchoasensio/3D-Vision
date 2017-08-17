function [ pointOut ] = computeTransformedPoints(H, point)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    pointOut = H * point;
    pointOut = pointOut / pointOut(3);
end

