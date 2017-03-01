function [ M ] = B_local( area,vertices )
%B_LOCAL Summary of this function goes here
%   Detailed explanation goes here
 r1 = vertices(1,1);
 r2 = vertices(2,1);
 r3 = vertices(3,1);
 M = area/(60)*[[6*r1+2*r2+2*r3 2*r1+2*r2+r3   2*r1+r2+2*r3];
                [2*r1+2*r2+r3   2*r1+6*r2+2*r3 r1+2*r2+2*r3];
                [2*r1+r2+2*r3   r1+2*r2+2*r3   2*r1+2*r2+6*r3]];
if(any(any(isnan(M))))
    disp('Nan during local B calculation')
end
end

