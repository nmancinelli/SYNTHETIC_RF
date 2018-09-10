function [SRF] = deconvolve(SRF, decon_method)
%Here we call Junlin's time-domain routine IDRF
%  below are the parameters he suggested to me in an email
%
    Parent = SRF.us;
    Daught = SRF.up;
    dt = SRF.time(2)-SRF.time(1);

    switch decon_method
        case 'IDRF'
            t_for = 5.;
            t_trun = 5.;
            gauss_t = 0.7;
            accept_mis = 0.00001;
            [RF_Time, RF] = DECONVOLUTION.IDRF(Parent,Daught,dt,t_for,t_trun,gauss_t,accept_mis);

        case 'ETMTM'
            TB=4;
            NT=3;
            winlen = 100;
            poverlap = 0.90;
            tag='data';

            %Locate max
            [~,imax] = max(abs(Parent));
            istart = round(imax-winlen/dt/2.);
            iend   = round(imax+winlen/dt/2.);
            Pwin = Parent(istart:iend)';
            Dwin = Daught(istart:iend)';

            [RF_Time, RF] = DECONVOLUTION.ETMTM(Pwin,Dwin,TB,NT,tag,dt,winlen,poverlap);
            [~, RFNORM]   = DECONVOLUTION.ETMTM(Pwin,Pwin,TB,NT,tag,dt,winlen,poverlap);
            RF = RF/max(abs(RFNORM));
    end
    
    %Flip polarity
    RF = -RF;
    
    SRF.rf_time = RF(RF_Time>=-100);
    SRF.time    = RF_Time(RF_Time>=-100);
    
end