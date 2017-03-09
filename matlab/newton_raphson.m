function [ u,v ] = newton_raphson( F,J,u0,v0,tol)
%NEWTON_RAPHSON Summary of this function goes here
%   Detailed explanation goes here
 load('vertices.dat');
 res = F(u0,v0);
 u = u0;
 v = v0;
 x = [u;v];
 N = size(u0,1);
figure
 while(norm(res,2)>tol)
     disp(norm(res,2))
     xmin = 0;
     xmax = 0.05;
     ymin = 0;
     ymax = 0.10;
     xlin = linspace(xmin,xmax,300);
    ylin = linspace(ymin,ymax,300);
    [X,Y] = meshgrid(xlin,ylin);
    U = griddata(vertices(:,1),vertices(:,2),u,X,Y,'linear');
    V = griddata(vertices(:,1),vertices(:,2),v,X,Y,'linear');

    subplot(1,2,1)
    mesh(X,Y,U)
    xlim([xmin xmax])
    ylim([ymin ymax])
    subplot(1,2,2)
    mesh(X,Y,V)
    xlim([xmin xmax])
    ylim([ymin ymax])
    drawnow;
     x = x-(J(u,v)\F(u,v));
     u = x(1:N);
     v = x(N+1:end);
     res = F(u,v);
 end
end

