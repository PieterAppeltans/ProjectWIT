function [ M ] = dRvdu( u,v,rq,Vmfv,Kmfu,Vmu,Kmu,Kmv )
%DRVDU Summary of this function goes here
%   Detailed explanation goes here
t= rq*dRudu(u,v, Vmu,Kmu,Kmv)+Vmfv*(-1./(Kmfu*(1+u/Kmfu).^2));
M =diag(t);

end

