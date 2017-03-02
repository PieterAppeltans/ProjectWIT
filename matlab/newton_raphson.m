function [ u,v ] = newton_raphson( F,J,u0,v0,tol)
%NEWTON_RAPHSON Summary of this function goes here
%   Detailed explanation goes here
 res = F(u0,v0);
 u = u0;
 v = v0;
 x = [u;v];
 N = size(u0,1);
 while(norm(res,inf)>tol)
     disp(norm(res,2))
     x = x-(J(u,v)\F(u,v));
     u = x(1:N);
     v = x(N+1:end);
     res = F(u,v);
 end
end

