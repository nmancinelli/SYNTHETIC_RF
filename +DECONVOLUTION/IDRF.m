function [RF_Time, RF] = IDRF(P,D,dt,t_for,t_trun,gauss_t,accept_mis)
% Iterative Deconvolution and Receiver-Function Estimation in time domain

%Hi,

%Here is the code I use.
%The t_for and t_trun are not important at all and I just make them to be 5.
%gauss_t is the pulse width which I use 0.7s. accept_mis I make it quite small,
%which controls when the iteration stops. The max iteration number is set
%in the function. P for parent, D for daughter and dt for sample rate.

%Junlin

t_max = t_for;

misfit = 1;
misfit_old = 99999;
misfit_ref = sqrt(sum(D.^2));

RF = zeros(length(P)*2-1,1);

D_cur = D;

itmax = 400;
%plot(D_cur)
%hold on
%for i=1:it_num

itnum = 0;
while misfit_old-misfit>accept_mis && itnum <= itmax
    [amp_corr,t_corr] = xcorr(D_cur,P);
    auto_corr = xcorr(P);
    [~,ind] = max(abs(amp_corr));
    amp_rf = amp_corr(ind)/auto_corr((length(t_corr)+1)/2);
    RF(ind) = RF(ind)+amp_rf;
    D_sub = conv(P,RF,'same');
    D_cur = D - D_sub;
    %plot(D_cur)
    misfit_old = misfit;
    misfit = sqrt(sum(D_cur.^2))/misfit_ref;
    itnum = itnum+1;
end


RF_Time = t_corr*dt;

RF(RF_Time>t_trun)=0;
RF = RF(RF_Time<=t_max);
RF_Time = RF_Time(RF_Time<=t_max);

if gauss_t~=0
    %gauss_sig = gauss_t/dt;
    % gauss_len = length(RF);
    % Gauss_win = gausswin(gauss_len,gauss_sig*4);
    gauss_sig = gauss_t/dt;
    x = linspace(-gauss_sig*4,gauss_sig*4,gauss_sig*8);
    Gauss_win = exp(-x.^ 2/(2*gauss_sig^2));
    
%     wl_T = 1/centfrq('mexh');
%     times_T = gauss_t/wl_T;
%     dt_wl = dt/times_T;
%     bd_wl = round(50/dt_wl)*dt_wl;
%     N_wl = 2*round(50/dt_wl)+1;
%     WL_win = mexihat(-bd_wl/2,bd_wl/2,N_wl);
%     n_fft_gau = 2^nextpow2(length(WL_win));
%     pad = ceil((n_fft_gau-length(WL_win))/2);
%     gau_tmp = [zeros(pad-1,1);WL_win';zeros(pad,1)];
%     gau_fft = fft(hilbert(gau_tmp));
%     gau_level = gau_fft(1);
%     gau_fft = fftshift(gau_fft);
%     gau_F_fft = 1/dt*((-n_fft_gau/2):(n_fft_gau/2-1))/n_fft_gau;
%     gau_fft = gau_fft./((1i*2*pi*gau_F_fft').^(1));
%     %gau_fft = -gau_fft./((1i*2*pi*gau_F_fft').^(3/2))./1i;
%     gau_fft = ifftshift(gau_fft);
%     gau_fft(1) = gau_level;
%     %RF_fft(1) = 0;
%     gau_shift = ifft(gau_fft);
%     if length(gau_shift) ~= length(WL_win)
%         gau_shift(end-pad+1:end) = [];
%         if pad>1
%             gau_shift(1:pad-1) = [];
%         end
%         
%         %gau_shift(length(Gauss_win)+1:end) = [];
%     end
%     gau_shift = real(gau_shift');
%     
%     RF = conv(RF,gau_shift,'same');
    RF = conv(RF,Gauss_win,'same');

end

RF = flipud(RF);
RF_Time = fliplr(RF_Time);

%plot(dt*t_corr(t_corr<0 & dt*t_corr>-40),RF_a(t_corr<0 & dt*t_corr>-40))

%pause



end
