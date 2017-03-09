function [ M ] = dRudv(u,v, Vmu,Kmu,Kmv)
%DRUDV Summary of this function goes here
%   Detailed explanation goes here
 t = (-(Vmu*u)./((Kmu+u).*(Kmv*(1+v/Kmv).^2)));
 M = diag(t);
end

