function [ M ] = A_local( area,vertices,Dr,Dz )
%A_LOCAL calculation of A for O2 (u) or CO2 (v) for one triangle

    G = [[(vertices(2,2)-vertices(3,2)) (vertices(3,1)-vertices(2,1))];
        [(vertices(3,2)-vertices(1,2)) (vertices(1,1)-vertices(3,1))];
        [(vertices(1,2)-vertices(2,2)) (vertices(2,1)-vertices(1,1))]];
    I = diag([Dr Dz]);
    M = (vertices(1,1)+vertices(2,1)+vertices(3,1))/(12*area)*G*I*G';

    if(any(any(isnan(M))))
        disp('Nan during local A calculation')
    end

end

