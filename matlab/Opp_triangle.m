function [ area ] = Opp_triangle( vertices)
%OPP_TRIANGLE Summary of this function goes here
%   Detailed explanation goes here
 area = abs(0.5*det([[vertices(2,1)-vertices(1,1) vertices(3,1)-vertices(1,1)];
                 [vertices(2,2)-vertices(1,2) vertices(3,2)-vertices(1,2)]]));
end

