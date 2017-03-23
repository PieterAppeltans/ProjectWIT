function [ M ] = dRvdv( u,v,rq,Vmu,Kmu,Kmv )
%DRVDV derivative of Rv to v

    M = rq*dRudv(u,v,Vmu,Kmu,Kmv);

end

