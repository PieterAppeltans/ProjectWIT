function [ M ] = C_local( length,vertices )
%C_LOCAL C for one mesh triangle

    r1 = vertices(1,1);
    r2 = vertices(2,1);
    M = length*[[r1/4+r2/12 r1/12+r2/12];[r1/12+r2/12 r1/12+r2/4]];
    
    if(any(any(isnan(M))))
        disp('Nan during local C calculation')
    end

end

