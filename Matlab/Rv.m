function [V] = Rv( u,v,rq,Vmfv,Kmfu,Vmu,Kmu,Kmv )
%RV Summary of this function goes here
%   Detailed explanation goes here
V = rq*Ru(u,v, Vmu,Kmu,Kmv)+Vmfv./(1+u./Kmfu);
if(isnan(any(V)))
    disp('Nan during Rv calculation')
end
end

