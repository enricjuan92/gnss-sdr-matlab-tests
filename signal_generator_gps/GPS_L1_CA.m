function [E1]=GPS_L1_Luis(Fs,BB_BW,CN0_IN,Tsig,FI,Doppler,Delay,numsat,flag_Datos,flag_local,flag_noise)
%[L1]=GPS_L1
%Fs= Sampling frequency [Hz]
%BB_BW= BaseBand bandwidth [Hz]
%CN0_IN= Carrier-to-noise-density ratio [dB]
%Tsig= Total signal time to generate [s]
%FI= Intermediate Frequency, 0 means baseband [Hz]
%Doppler = Doppler frequency [Hz]
%Delay = Code delay [s]
%numsat= Satellite Vehicle number(1-31)
%flag_datos = 1 -> Signal contains random telemetry bits 0 -> No telemetry
%flag_local = 1 -> Signal will be used as a local replica
%flag_noise = 1 -> Signal contains noise 0 -> signal is noise-free
%%*** generador se�al BPSK GPS L1 C/A***
%***** Javier Arribas L�zaro CTTC 2010 ***
%**** v 1.0 ****
%% Definiciones:
BB_SAT=20E6;
%Chip Rate
Rc_L1_CA=1.023e6; %Rc,E1-B

%Symbol Rate (NAV MESSAGES)
Rd_L1_DATA=50;

%% Generacion se�ales PRN

num_chips=floor(Tsig/(1/Rc_L1_CA)); %numero de chips para L1 CA en Tsig
C_L1=digitGPS_L1_CA(num_chips,Rc_L1_CA,0,numsat,1); %generaci�n de se�al PRN C_E1_B muestreada a Rc_E1_B (directamente chip rate)

%**** chip DATOS generation (FAKE) 
numdatos=ceil(Tsig/(1/Rd_L1_DATA));
datos_nav=round(rand([numdatos,1])); %secuencia de datos para la entrada I, pero hay que extenderla!
datos_nav=datos_nav.*-2.+1;

if numdatos==1
    D_L1=resample(datos_nav,Rc_L1_CA,Rd_L1_DATA,0); %usaremos resample para hacer un "upsample"
else
    D_L1=resample(datos_nav,Rc_L1_CA,Rd_L1_DATA,0)'; %usaremos resample para hacer un "upsample"
end

%delay_samples=floor(Delay*Fs);
%D_L1=cshift(D_L1,delay_samples);

if flag_Datos==1
    nsamples=length(C_L1);
    e_L1=C_L1.*D_L1(1:nsamples); %c�lculo de la se�al C_E1_B*D_E1_B a Fs = 1/chip_rate   PRN(I)
else
    e_L1=C_L1;
end

n_samples=floor(Tsig/(1/Fs));
t=(1:n_samples).*(1/Fs);

S_E1=sampledCode(e_L1,Delay,t,1/Rc_L1_CA,1/Fs);

if flag_local==0 % Satellite signal generation
    n_samples=floor(Tsig/(1/Fs));
    numSamplesPsig=floor(1e-3/(1/Fs)); 
    
    % filtrado para BB_SAT
    if BB_BW>BB_SAT
        [b,a] = butter(5,BB_SAT/(Fs/2),'low');
        S_E1_filt=filtfilt(b,a,S_E1); %filtfilt tiene retardo cero pero eleva
    else
        [b,a] = butter(5,BB_BW/(Fs/2),'low');
        S_E1_filt=filtfilt(b,a,S_E1); %filtfilt tiene retardo cero pero eleva
    end
    
    %al cuadrado la respuesta del filtro..
    Pclean=sum(abs(S_E1_filt(1,1:numSamplesPsig)).^2)/numSamplesPsig;
    P_signal=(10^(CN0_IN/10))*(1/(BB_BW*2));
    S_E1_BB=(1/sqrt(Pclean)).*sqrt(2*P_signal).*S_E1_filt; %se�al normalizada y con la amplitud adecuada para la C/N0 considerando Pn=1 (ruido blanco normalizado)
    
    if flag_noise==1
        %% *** AGW Noise ***
        % consideramos el BW de ruido en pasa banda, por lo tanto la densidad
        % de ruido se considerar� tomando BW=2xBB_BW
        ruido=randn(1,length(S_E1_BB))+randn(1,length(S_E1_BB)).*i;
        [b,a] = butter(5,BB_BW/(Fs/2),'low');
        ruido_filt=filtfilt(b,a,ruido);
        p_ruido_filt=(ruido_filt(1:numSamplesPsig)*ruido_filt(1:numSamplesPsig)')/numSamplesPsig;
        ruido_filt=(sqrt(2)/sqrt(p_ruido_filt))*ruido_filt;
        E1_BB=S_E1_BB+ruido_filt; %a�adimos ruido filtrado a FI
        
        % DEBUG: SNR and CN0 MEASUREMENT
            Ps_lin=(S_E1_BB(1,1:numSamplesPsig)*S_E1_BB(1,1:numSamplesPsig)')/numSamplesPsig
            Pn_lin=(ruido_filt(1,1:numSamplesPsig)*ruido_filt(1,1:numSamplesPsig)')/numSamplesPsig
            SNR=Ps_lin/Pn_lin % Carrier to Noise ratio
            SNR_db=10*log10(SNR) % Carrier to Noise ratio (dB)
            CN0=10*log10(SNR*(BB_BW*2))
        
    else
        E1_BB=S_E1_BB;
    end

    %% ****** La subimos a FI *****
    %**** L.O signal generation ****
    LO_freq=FI+Doppler;
    if LO_freq~=0 %mod 13/9/2010, antes LO_freq>0
        t=(1:n_samples).*(1/Fs);

        l_osc=exp(2*pi*LO_freq*t*j); % subimos a FI %mod 13/9/2010, now is baseband, antes l_osc=sin(2*pi*LO_freq*t);

        E1=E1_BB.*l_osc; %mezclamos
    else
        E1=E1_BB;
    end
else %local code generator
    E1=S_E1; %clean signal, baseband without filter
end