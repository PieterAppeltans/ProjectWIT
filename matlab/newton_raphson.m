function [ u,v ] = newton_raphson( F,J,u0,v0,tol )
%NEWTON_RAPHSON Newton-Raphson for system of nonlinear equations

    load('vertices.dat');
    res = tol+1;
    u = u0;
    v = v0;
    x = [u;v];
    N = size(u0,1);
    it = 0;
    while(res>tol && it<150)
        disp(res)
        x_prev = x;
        x = x-(J(u,v)\F(u,v));
        u = x(1:N);
        v = x(N+1:end);
        res = norm(x-x_prev);
        it = it+1;
    end

end
