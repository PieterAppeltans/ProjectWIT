load('boundary.dat');
load('vertices.dat');
figure
hold on
for i=1:size(boundary,1)
    v1 =vertices(boundary(i,1),:);
    v2 = vertices(boundary(i,2),:);
    plot([v1(1),v2(1)],[v1(2),v2(2)],'b-');
end