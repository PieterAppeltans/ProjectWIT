function [ M ] = dRudv(u,v, Vmu,Kmu,Kmv)
%DRUDV Summary of this function goes here
%   Detailed explanation goes here
 t = ((Vmu*u)./(Kmu*u)).*(-1./(Kmv*(1+v/Kmv)));
 M = diag(t);
end

