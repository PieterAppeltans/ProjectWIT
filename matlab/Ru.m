function [ V ] = Ru( u,v, Vmu,Kmu,Kmv )
%RU Calculate Ru value

    V = (Vmu*u)./((Kmu+u).*(1+v/Kmv));
    if(any(isnan(V)))
        disp('Nan during Ru calculation')
    end

end

