function [ M ] = dRudu( u,v, Vmu,Kmu,Kmv )
%DRUDU Summary of this function goes here
%   Detailed explanation goes here
 t = Vmu*Kmu./((1+v./Kmv)*(Kmu+u).^2);
 diag(t);
end

