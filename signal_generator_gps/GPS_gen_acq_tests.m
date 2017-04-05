% /*!
%  * \file GPS_gen_acq_tests.m
%  * \brief GPS L1/CA signal generator for several satellites
%
%
%  * \author Luis Esteve, 2012. luis(at)epsilon-formacion.com
%  * -------------------------------------------------------------------------

clc;
clear all;
close all;

%% Satellite signals parameters

Fs= 4e6; %Sample frequency (Hz)
BB_BW=Fs/2.1; % Bandwidth (Hz) (the parameter has to be slightly below Fs/2 
              % to avoid instability problems)
CN0_Sats= 48; %CN0 for all the satellites
Tsig=200e-3; %Time interval to be generated (s)

FI=0;    %Intermediate Frequency (Hz)
n_samples=Fs*Tsig;

flag_Datos = 0; % Data flag (0=no data bits, 1=data bits simulated)
flag_local = 0; % Local/transmitted signal (0=local signal, not filtered
                % 1 = transmitted signal from satellite, filtered)
flag_noise = 1; % Adding channel noise (0=without noise, 1=with noise)
N_sats = 2; % Number of signal satellite to be generated
ID_sat = [1 11 14 20 22 23 31 32]; % ID's of the satellites
Doppler = [3952 4906 2792 8163 2663 8125 5501 3868]; % Doppler of each signal (Hz)
Delay_sps = [3767 3037 3503 1865 3698 1062 2407 79]; % Code delay of each signal (samples)
Delay_sec = Delay_sps/Fs;
m_seconds = 20000;
SIGMA2_NOISE=1; %this value should be 1 to keep satellites CN0...

% Generation of the signal of the first satellite

s=GPS_L1_CA(Fs,BB_BW,CN0_Sats,Tsig,FI,Doppler(1),Delay_sec(1),ID_sat(1),flag_Datos,flag_local,flag_noise);




if (N_sats > 1)
    
% Generation of the signals to complete the set of signals to be generated:
% flag_noise is set to zero because the noise was added in the previous signal
% 
    
    flag_noise = 0;

    for k=2:1:N_sats
               s=s+GPS_L1_CA(Fs,BB_BW,CN0_Sats,Tsig,FI,Doppler(k),Delay_sec(k),ID_sat(k),flag_Datos,flag_local,flag_noise);
    end
end


name = strcat('/media/DATA/Proyectos/GSoC/2015/signals/GPS_L1_CA_',num2str(Fs),'_sps_CN0_', num2str(CN0_Sats),'.dat');

write_apend_complex_binary(s, name, n_samples);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternative way to generate the band pass noise
% Only use with gflag = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P_signal=(10^(CN0_Sats/10))*(1/(BB_BW*2))
% Ps = (s*s')/n_samples
% 
% t=(1:n_samples).*(1/Fs);
% l_osc=exp(2*pi*FI*t*1i);
% [b,a] = butter(5,BB_BW/(Fs/2),'low'); %Butterworth Low Pass Filter
% 
% for k=1:1:(m_seconds/20)
%     
%     noise=sqrt(SIGMA2_NOISE)*(randn(1,n_samples)+randn(1,n_samples)*1i);
%     noise_filt=filtfilt(b,a,noise);
%     p_noise_filt=(noise_filt(1:n_samples)*noise_filt(1:n_samples)')/n_samples;
%     noise_filt=(sqrt(2)/sqrt(p_noise_filt))*noise_filt;
%     band_pass_noise = noise_filt.*l_osc;
% 
%     Pn = ( band_pass_noise* band_pass_noise')/n_samples
% 
%     SNR_lin = Ps/Pn
%     SNR_db=10*log10(SNR_lin) 
%     
% % Carrier to Noise ratio (dB)
%     CN0=10*log10(SNR_lin*(BB_BW*2))
%     x = s +  band_pass_noise;
% 
%     write_apend_complex_binary(x, name, n_samples);
%     
%     clear noise;
%     clear noise_filt;
%     clear noise_pass_band;
%     clear x;
% end



 