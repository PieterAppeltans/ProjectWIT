function [ M ] = dRvdv( u,v,rq,Vmu,Kmu,Kmv )
%DRVDV Summary of this function goes here
%   Detailed explanation goes here

M = rq*dRudv(u,v,Vmu,Kmu,Kmv);

end

