function [ V ] = Ru(u,v, Vmu,Kmu,Kmv)
%RU Summary of this function goes here
%   Detailed explanation goes here
V = (Vmu*u)./((Kmu+u).*(1+v/Kmv));
if(isnan(any(V)))
    disp('Nan during Ru calculation')
end
end

