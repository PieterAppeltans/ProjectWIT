function [V] = Rv( u,v,rq,Vmfv,Kmfu,Vmu,Kmu,Kmv )
%RV Calculate value of Rv

    V = rq*Ru(u,v, Vmu,Kmu,Kmv)+Vmfv./(1+u./Kmfu);
    if(any(isnan(V)))
        disp('Nan during Rv calculation')
    end
    
end

