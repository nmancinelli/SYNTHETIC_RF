function [SRF] = compute_synthetic_wvfrms(SRF)
    zt = SRF.VelocityModel.zlayt;
    zb = SRF.VelocityModel.zlayb;
    vp = SRF.VelocityModel.vplay;
    vs = SRF.VelocityModel.vslay;
    rh = SRF.VelocityModel.rhlay;
    RayParam = SRF.rayparam;
    PERIOD = SRF.period;

    [SRF.ux, ~, SRF.uz, SRF.time] = PROPMAT.write_propmat_syn(zt,zb,vp,vs,rh,RayParam,PERIOD);

end