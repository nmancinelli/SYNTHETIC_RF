function [SRF] = deconvolve(SRF)
%Here we call Junlin's time-domain routine IDRF
%  below are the parameters he suggested to me in an email
%
    t_for = 5.;
    t_trun = 5.;
    gauss_t = 0.7;
    accept_mis = 0.00001;
    Parent = SRF.us;
    Daught = SRF.up;
    dt = SRF.time(2)-SRF.time(1);

    [RF_Time, RF] = DECONVOLUTION.IDRF(Parent,Daught,dt,t_for,t_trun,gauss_t,accept_mis);
    
    %Flip polarity
    RF = -RF;
    
    SRF.rf_time = RF(RF_Time>=-100);
    SRF.time    = RF_Time(RF_Time>=-100);
    
end