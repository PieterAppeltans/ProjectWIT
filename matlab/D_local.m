function [ V ] = D_local( length,vertices)
%D Summary of this function goes here
%   Detailed explanation goes here

r1 = vertices(1,1);
r2 = vertices(2,1);
V = length*[r1/3.+r2/6.;r1/6+r2/3];
if(any(isnan(V)))
    disp('Nan during local D calculation')
end
end

