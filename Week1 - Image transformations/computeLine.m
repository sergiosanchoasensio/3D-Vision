function [ line ] = computeLine( pointOne, pointTwo)
%computeLine Summary of this function goes here
%   Detailed explanation goes here
a = pointOne(2) - pointTwo(2);
b = pointTwo(1) - pointOne(1);
c = pointOne(1)*pointTwo(2) - pointTwo(1)*pointOne(2); 

line = [a, b, c];
end

