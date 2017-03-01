function [ l ] = Length_edge( vertices )
%LENGTH_EDGE Summary of this function goes here
%   Detailed explanation goes here
 l =norm(vertices(1,:)-vertices(2,:),2);

end

