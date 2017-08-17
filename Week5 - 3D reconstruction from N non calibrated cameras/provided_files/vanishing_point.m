function [v] = vanishing_point(xo1, xf1, xo2, xf2)
% computes the vanishing point formed by the line 
% that joins points xo1 and xf1 and the line 
% that joins points x02 and xf2
    v = cross(cross(xo1, xf1), cross(xo2, xf2));
end
