function [ M ] = dRudu( u,v,Vmu,Kmu,Kmv )
%DRUDU derivative of Ru function to u

    M = diag((Vmu*Kmu)./((1+v/Kmv).*(Kmu+u).^2));
 
end

