function [ u,v ] = newton_raphson( F,J,u0,v0,tol)
%NEWTON_RAPHSON Summary of this function goes here
%   Detailed explanation goes here
 load('vertices.dat');
 res = tol+1;
 u = u0;
 v = v0;
 x = [u;v];
 N = size(u0,1);
 it = 0;
 while(res>tol && it<50)
     disp(res)
     x_prev = x;
     x = x-(J(u,v)\F(u,v));
     u = x(1:N);
     v = x(N+1:end);
     res = norm(x-x_prev)/norm(x);
     it = it+1;
 end
end

