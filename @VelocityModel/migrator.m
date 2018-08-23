function [Z,T] = migrator(VelocityModel,ray_parameter)
    
    nlay = length(VelocityModel.zlayb);
    masked = false;

    function [TInc] = timeIncrement(vp,vs,ztop,zbot)
        A = sqrt(vs^-2 - ray_parameter^2);
        B = sqrt(vp^-2 - ray_parameter^2);
        h = zbot - ztop;
        TInc = (A-B)*h;
        if (vp^-2 - ray_parameter^2) < 0 && masked == false;
            masked = true;
        end
    end

    t = 0;
    T=zeros(nlay,1);
    for ilay = 1:nlay;
        [tinc] = timeIncrement(VelocityModel.vplay(ilay),VelocityModel.vslay(ilay),VelocityModel.zlayt(ilay),VelocityModel.zlayb(ilay));
        
        if masked;
            T(ilay) = NaN;
        else
            t = t + tinc;
            T(ilay) = t;
        end

    end
    
    Z = VelocityModel.zlayb;
    T = [0;T];
    Z = [0;Z];
        

end