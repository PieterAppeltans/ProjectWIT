function [ M ] = A_local( area, vertices,Dr,Dz)
%A_LOCAL Summary of this function goes here
%   Detailed explanation goes here
G = [[sqrt(Dr)*(vertices(2,2)-vertices(3,2)),sqrt(Dz)*(vertices(3,1)-vertices(2,1))];
    [sqrt(Dr)*(vertices(3,2)-vertices(1,2)),sqrt(Dz)*(vertices(1,1)-vertices(3,1))];
    [sqrt(Dr)*(vertices(1,2)-vertices(2,2)),sqrt(Dz)*(vertices(2,1)-vertices(1,1))]];
M = (vertices(1,1)+vertices(2,1)+vertices(3,1))/(12*area)*G*G';
if(any(any(isnan(M))))
    disp('Nan during local A calculation')
end
end

