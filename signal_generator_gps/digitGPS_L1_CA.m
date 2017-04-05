function code2 = digitGPS_L1_CA(n,fs,offset,numsat,level)
%[code2]=digitCE1B[n,fs,offset,numsat,level]
%Genera la señal PRN C_E1B del satélite Galileo especificado
%n es el numero de chips a generar
%fs es la frecuencia de muestreo
%offset es el desplazamiento temporal inicial
%numsat es el número de satélite
%level en caso de ser =1 los niveles de salida serán +-1
%en caso contrario serán +1 y 0
%****JAVIER ARRIBAS LÁZARO **********

symb_rate = 1.023e6; %Symbol rate en Hz
%*** primero vamos a cargar el código PRN del sat
code=codegen(numsat);

%*** empezamos a copiar con offset ***
code_ext=zeros(1,n);
l=1 + offset; %empezamos con el offset
for k=1:1:n
   code_ext(k)=code(l);
   l=l+1;
   if l>1023
       l=1;
   end
end

if fs>symb_rate
    %*** upsampling a Fs ****
    prn=resample(code_ext,fs,symb_rate,0); %usaremos resample para hacer un "upsample"
else
    prn=code_ext;
end
%*** cambiamos de niveles si se solicita **

if level==0
    prn=(prn+1)./2; %pasamos a niveles 0 +1 (1 -> -1 0 -> +1
end

code2=prn;
end

