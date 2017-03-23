function [ M ] = dRudv( u,v,Vmu,Kmu,Kmv )
%DRUDV derivative of Ru to v

    M = diag(-(Vmu*u)./((Kmu+u).*(Kmv*(1+v/Kmv).^2)));
 
end

