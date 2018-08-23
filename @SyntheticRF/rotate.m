function [SRF] = rotate(SRF)
    a = SRF.VelocityModel.vplay(1);
    b = SRF.VelocityModel.vslay(1);
    rp = SRF.rayparam;
    
    qa = sqrt((a^-2)-(rp^2) );
    qb = sqrt((b^-2)-(rp^2) );
    A  = rp*(b^2)/a;
    B  = ((b^2)*(rp^2)-0.5)/(a*qa);
    C  = (0.5-(b^2)*(rp^2))/(b*qb); 
    D  = rp*b;
    
    %flip z polarity
    SRF.up  =  A*SRF.ux + B*-SRF.uz;
    SRF.us  =  C*SRF.ux + D*-SRF.uz;
    
end
