function [ M ] = A_local( area, vertices,Dr,Dz)
%A_LOCAL Summary of this function goes here
%   Detailed explanation goes here
G = [[(vertices(2,2)-vertices(3,2)) (vertices(3,1)-vertices(2,1))];
    [(vertices(3,2)-vertices(1,2)) (vertices(1,1)-vertices(3,1))];
    [(vertices(1,2)-vertices(2,2)) (vertices(2,1)-vertices(1,1))]];
I = diag([Dr Dz]);
M = (vertices(1,1)+vertices(2,1)+vertices(3,1))/(12*area)*G*I*G';
if(any(any(isnan(M))))
    disp('Nan during local A calculation')
end
end

