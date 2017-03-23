function [ M ] = dRvdu( u,v,rq,Vmfv,Kmfu,Vmu,Kmu,Kmv )
%DRVDU derivative of RV to u

    M = rq*dRudu(u,v, Vmu,Kmu,Kmv)+diag(Vmfv*(-1./(Kmfu*(1+u/Kmfu).^2)));

end

