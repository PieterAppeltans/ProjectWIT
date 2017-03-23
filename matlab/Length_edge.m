function [ l ] = Length_edge( vertices )
%LENGTH_EDGE calculates length of edge

    l =norm(vertices(1,:)-vertices(2,:),2);

end

