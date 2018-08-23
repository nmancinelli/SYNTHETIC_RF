function [SRF] = migrate_to_depth(SRF)
    depths = 0:300;
    npts = length(depths);
    
    [Z,T] = SRF.VelocityModel.migrator(SRF.rayparam);
    
    for ipt = 1:npts;
        depth = depths(ipt);
        time = -interp1(Z,T,depth);
        if isnan(time)
            value = nan;
        else
            value = interp1(SRF.time, SRF.rf_time, time);
        end
        SRF.rf_depth(ipt) = value;
        %fprintf('%f %f %f !!!!\n',depth,time,value)
        
    end
    SRF.depth = depths;
end