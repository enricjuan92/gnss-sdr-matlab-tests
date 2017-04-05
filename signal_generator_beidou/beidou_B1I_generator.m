% /*!
%  * \file beidou_B1I-generator.m

%  create a file containing the desired B1I code sequence of the selected satellite 
%
%
%  * -------------------------------------------------------------------------

clc;
clear all;
close all;


%% ------------Settings---------------------------------------------------------
fileName1 = 'FFF005.dat';
% fileType = 1;
fileType = 2;

f_rf = 1561.098e6;    %[Hz] BeiDou B1 nominal frequency;
f_if = 0.098e6;       %[Hz] IF nominal frequency;
f_prn= 2.046e6;       %[Hz] Nominal PRN-generator clock frequency;
f_nh = 1e3;           %[Hz] Nominal Neiman-Huffman-generator clock frequency;
f_data = 50;          %[Hz] Nominal data rate of BeiDou D1 signal
f_data_NH = 1e3;      %[Hz] Data bit rate after modulation by NH code

k_car_prn  = f_rf / f_prn;            %[unitless] Ratio between RF frequency and PRN clock freq;
k_car_nh   = f_rf / f_nh;             %[unitless] Ratio between RF frequency and NH clock freq;
k_car_data = f_rf / f_data;           %[unitless] Ratio between RF frequency and data clock freq;
k_car_data_NH = f_rf / f_data_NH;     %[unitless] Ratio between RF frequency and data clock freq after modulation;

phi0_if   = 0;      %[rad] Initial phase of RF signal;
phi0_prn1 = 0;      %[rad] Initial phase of PRN signal;

phi0_nh   = 0;      %[rad] Initial phase of NH signal;
phi0_bk   = 0;      %[rad] Initial phase of BK signal;
phi0_data = 0;      %[rad] Initial phase of data signal;

f_d = 2800;         %[Hz] Initial Doppler frequeny for RF-signal;
df  = -0.55;        %[Hz/sec] Initail Doppler frequency change rate for RF-signal; ASK

fs = 16.00e6;      %[Hz] Sampling frequency;
ts = 1/fs;         %[sec]

T  = 4;            %[sec] Signal length to be generated;
% T = 200e-3;      %  This is the Time interval fot he signal_generator for GPS

T_elem = 10e-3;    %[sec] The smallest signal part to be generated.
T_parts = T/T_elem;%[unitless] On how many segments T will be divided.
dT = T / T_parts;  %[sec]

prn_num  = 7;        %PRN number;
prn_len  = 2046;     %[chips] PRN-code length in chips (bits);
nh_len   = 20;       %Heiman-Haffman code length in chips (bits);
data_len = 36000;    %Navigation message length in bits (after convolutional coder); 


%Noise generation is under development...
%SNR_dB = -3;             %SNR of the generating signal;
%SNR = 10^(SNR_dB/10);    %SNR in times;
%noise_pwr = 1/SNR;       %noise power for predefined SNR;
%noise_pwr = noise_pwr * ( (fs/2) / (2*f_prn) ); %Take into accound sampling frequency and signal bandwidth!



%% -----------Signal_generator--------------------------------------------------

t = ts : ts : T_elem; %time samples for generating signal of T_elem length;

PRN1 = generateB1Icode(prn_num);                                      %Generate PRN-code;
NH_original = [0 0 0 0 0 1 0 0 1 1 0 1 0 1  0  0  1  1  1  0];        %Generate Neiman-Huffman-code
NH = [-1 -1 -1 -1 -1 1 -1 -1 1 1 -1 1 -1 1 -1 -1  1  1  1  -1];
DATA = (2 * round(rand(1, data_len))) - 1;                            % Temporary Stub for data bits after convolutional coder; (pseudorandom values of 1 and -1, of data_len values)

signal_I = [];
%signal_Q = [];

[fd1, err1] = fopen(fileName1, 'ab');

for k=1:T_parts
    
    %Calculate phase for carrier;
    phi_if = (phi0_if) + (2*pi*f_if*t) + (2*pi*f_d*t);               % INITIAL PHASE + IF FREQ + DOPPLER FREQ

    %phi_if = (phi0_if) + (2*pi*f_if*t) + (2*pi*f_d*t) + ((2*pi*df*t).*t);
    %phi0_if            - this is initial phase;
    %2*%pi*f_if*t       - this is phase change due to nominal frequency;
    %2*%pi*f_d*t        - this is phase change due to Doppler (f_d) frequency;
    %(2*%pi*df*t).*t    - this is phase change due to Change of Doppler (acceleration of the object);
    %t must be replaced by (t + (k-1)*T_elem)! (Must it? ;) )


    phi_if = mod(phi_if, (2*pi));                                   % Convert carrier phase to the range [0..2*pi];  
    phi0_if = phi_if(end);

    %f_d = f_d + df*T_elem;                                         % FOR THE MOMENT I CONSIDER THE DOPPLER CONSTANT IN TIME (NO ACCELLERATION OF THE OBJECT)
    
    %Calculate phase for PRN1;
    phi_prn1 = (phi0_prn1) + ((f_rf/k_car_prn)*t) + ( (f_d/k_car_prn)*t );
    %phi_prn1 = (phi0_prn1) + ((f_rf/k_car_prn)*t) + ( (f_d/k_car_prn)*t ) + ( ((f_d/100)/k_car_prn)*t ) + ((( (df/k_car_prn) *t).*t));
    %The change in phase exactly the same like for carrier. The only difference is a special multiplier "k_car_prn".
%     phi0_prn1 = modulo( phi_prn1(end), prn_len );   
    
    phi0_prn1 = mod( phi_prn1(end), prn_len );         % from modulo (scilab) to mod (matlab)
    prn1_indx = (fix(mod(phi_prn1, prn_len)));
    
    %Calculte phase for Neiman-Huffman;
    phi_nh =  (phi0_nh) + ((f_rf/k_car_nh)*t) + ( (f_d/k_car_nh)*t );

    %phi_nh =  (phi0_nh) + ((f_rf/k_car_nh)*t) + ( (f_d/k_car_nh)*t ) + ((( (df/k_car_nh) *t).*t));
    phi0_nh = mod(phi_nh(end), nh_len );
    nh_indx = (fix( mod(phi_nh, nh_len) ));
     
    %Calculate phase for data-message;
    phi_data =  (phi0_data) + ((f_rf/k_car_data)*t) + ((f_d/k_car_data)*t) + ((( (df/k_car_data) *t).*t));
    phi0_data = mod(phi_data(end), data_len);
    data_indx = (fix(mod(phi_data, data_len)));
    
    %Generate carrier;
%     carr_sin = sin(phi_if);    %generate carrier (I)
%     carr_cos = cos(phi_if);    %generate carrier (Q)    % FOR GLONASS

%     carr_sin = sin(phi_if);    %generate carrier (Q)
    carr_cos = cos(phi_if);    %generate carrier (I)    % FOR BEIDOU

    %Generate PRN;
    prn1 = PRN1(prn1_indx+1);    %generate PRN for I-channel;
%     prn2 = PRN2(prn2_indx+1);    %generate PRN for Q-channel;
    
    %Generate NH;
    nh = NH(nh_indx+1);         %generate NH-code for pilot-channel;
    
    
    %Generate DATA;
    data = DATA(data_indx+1);
    
    if (fileType == 1) %TO DO: make this it work (fileType=1)!
      s_I =   ( (prn1 .* nh).* data  ) ;                         % I moduled the B1I with the data
      %s_Q = (  (prn2 .* bk) .* data  );
      signal_RSLT = (1/sqrt(2))*(s_I.*carr_sin) + (1/sqrt(2))*(s_Q.*carr_cos);
      
    else
      s_I =    ( (prn1 .* nh).* data  ) ;
      %s_Q = (  (prn2 .* bk) .* data  );
      
      %The following commented lines are left in order to understand how the signal is generated!
      %signal_RSLT(1:2:2*length(s_I)-1) =   ((s_I .* carr_sin) - (s_I .* carr_cos));
      %signal_RSLT(2:2:2*length(s_Q))   =   ((s_I .* carr_sin) + (s_I .* carr_cos));
      %signal_RSLT(1:2:2*length(s_I)-1) =   ((s_Q .* carr_sin) + (s_Q .* carr_cos));
      %signal_RSLT(2:2:2*length(s_Q))   =   -((s_Q .* carr_sin) - (s_Q .* carr_cos));
      %signal_RSLT(1:2:2*length(s_I)-1) =   ((s_I .* carr_sin) - (s_I .* carr_cos)) + ((s_Q .* carr_sin) + (s_Q .* carr_cos));
      %signal_RSLT(2:2:2*length(s_Q))   =   ((s_I .* carr_sin) + (s_I .* carr_cos)) - ((s_Q .* carr_sin) - (s_Q .* carr_cos));
      
      %After simplification the result formula looks like:
      signal_RSLT(1:2:2*length(s_I)-1) = ((s_I ).*carr_cos);
%       signal_RSLT(2:2:2*length(s_Q))   = ( (s_I - s_Q).*carr_sin ) + ( (s_I + s_Q).*carr_cos );
      
      %Next two lines are wrong! They are left as an example of
      %how not to do!!!
      %signal_RSLT(1:2:2*length(signal_I)-1) =  signal_I;
      %signal_RSLT(2:2:2*length(signal_Q))   = -signal_Q;
    end
    
    %Let's add some noise:
    %Noise generation is under development...
    %noise = grand(1,length(signal_RSLT),"nor",0,noise_pwr);
    %signal_RSLT = signal_RSLT + noise;
    %signal_RSLT = signal_RSLT / max(abs(signal_RSLT)); %Normalization after adding noise;
    
    %Next step is to pass signal througn bandpass filter;
    %TODO: make bandpass filter;
    
    %Next step is to write signal to file;
    %Can't understand the fact that when using 256-level discretization - problems uccure!
    %With 4-level discretization everything works fine!
    %signal_RSLT_r = (   round( (63.5 * signal_RSLT) + 63.5) * 2   ) - 127; //Discretization of quasianalogue signal;
    signal_RSLT_r = (round( (1.5 * signal_RSLT) + 1.5) * 2) - 3; %Discretization of quasianalogue signal;
    v=fwrite(fd1, signal_RSLT_r, 'int');
    
end

fclose(fd1);

