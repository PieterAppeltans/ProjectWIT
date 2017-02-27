function [ M ] = dRvdv( u,v,rq,Vmfv,Kmfu,Vmu,Kmu,Kmv  )
%DRVDV Summary of this function goes here
%   Detailed explanation goes here
t = rq*dRudv(u,v, Vmu,Kmu,Kmv);
M = diag(t);

end

