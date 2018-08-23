classdef SyntheticRF
    %SyntheticRF A synthetic receiver function
    %   Codes to create one from a velocity model
    
    properties
        rayparam
        period
        VelocityModel
        rf_time
        time
        rf_depth
        depth
        masking_depth
        uz
        ux
        us
        up
    end
    
    methods
        function SRF = SyntheticRF(VelocityModel,ray_parameter, period)
            SRF.VelocityModel = VelocityModel;
            SRF.rayparam = ray_parameter;
            SRF.period = period;
        end
    end
    
end