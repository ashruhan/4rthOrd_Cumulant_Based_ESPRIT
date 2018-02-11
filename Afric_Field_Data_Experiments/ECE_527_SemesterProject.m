function ECE_527_SemesterProject

global handles 

figure(handles.window)

NUM_Targets = (str2num(get(handles.NUM_Targets,'String')));
Band_WDT = (str2num(get(handles.Band_WDT,'String')));
Trans_Power = (str2num(get(handles.Trans_Power,'String')));
Pulse_Width = (str2num(get(handles.Pulse_Width,'String')));
Noise_Pwr = (str2num(get(handles.Noise_Pwr,'String')));
NOM_Range = (str2num(get(handles.NOM_Range,'String')));
Ant_Gain = (str2num(get(handles.Ant_Gain,'String'))); %Number of tagets  
Smearing = (str2num(get(handles.Smearing,'String')));
% Global declarations that the user will input 
% keeping at a constant untill matched filter works and plots correctly

      
     
Band_WDT=Band_WDT*10^6;  %Chirp Bandwidth
Trans_Power=Trans_Power*10^6;%Transmitted Power
Pulse_Width=Pulse_Width*10^-6; %Pulse width MUST SATISFY NYQUIST RATE 
Noise_Pwr=Noise_Pwr*10^-3 ; %Noise Power
NOM_Range=NOM_Range*10^3;     %Nominal Rang
Ant_Gain=exp(Ant_Gain)/10;           %Antenna Gain User input will be in dB so account for that



% Range Parameters
thous=10000;

c=300*10^6;

f=1000*10^6;%frequency

lambda=c/f;%wavelength

beta=(Band_WDT)/(Pulse_Width*2);% Beta value for S

betaPrimed = Smearing*(Band_WDT)/(Pulse_Width*2);% beta values for H

dt=1/(2*Band_WDT); % Time Domain Sampling Interval

R = NOM_Range + 9000*rand(1,NUM_Targets);%Random targets Range of the platform relys on N  
%R=[8500,2000];

N = round(Pulse_Width/dt); NN = N/2;

t = ([1:N]-NN)*dt; %time of sampling

tn=((2.*R)./c);%time shift



%%%%%%%%%%%%%%%%%%%%%%%Return Echo summation%%%%%%%%%%%%%%%%

RCS=1+9*rand(1,NUM_Targets); %Random Radar cross section relys on N

Pr=((Trans_Power.*Ant_Gain.^2.*lambda.^2.*RCS)./(((4*pi).^3).*R.^4)); %Radar Equation

Rn=sqrt(Pr).*exp(1i.*2.*pi.*rand(1,NUM_Targets));

G = exp(-1i.*2.*pi.*beta.*t.^2);%Linear Chirp Equation

H=exp(1i.*2.*pi.*betaPrimed.*t.^2);%Matched Filter

   Rx_echo = sqrt(-Noise_Pwr*sqrt(1-rand(1,thous))).*exp(1i*2*pi*rand(1,thous));
   
for i=1:NUM_Targets  
    
    M = round(tn(i)/(.65.*(Band_WDT/(50*10^6))*dt));
    
    Rx_echo(M-NN:M+NN-1) = Rn(i).*G(1:N);
   
end

F=filter(H,NUM_Targets,Rx_echo); 

%%%%%%%%%%%%%%%%%%%%%%%Return Echo summation%%%%%%%%%%%%%%%%

Range=zeros(1,length(F));
for j=1:NUM_Targets;
    Range(round(R(j)))=1;
end
time=zeros(1,length(F));
for j=1:length(time);
   time(j)=(1+Band_WDT/10^8)*j; 
end
subplot(3,1,1)
plot(Range)
xlim([NOM_Range NOM_Range+9000])
axes('position', [30,30, 60,60])
xlabel('Range');
ylabel('Unitless');
subplot(3,1,2)

plot(real(Rx_echo))
xlim([NOM_Range NOM_Range+9000])
xlabel('Range');
ylabel('Power Recieved')

subplot(3,1,3)
plot(abs(F))
xlim([NOM_Range NOM_Range+9000])
xlabel('Range');
ylabel('Image Amplitude')

subplot(3,1,4)

cla(findall(handles.window,'type','axes'));

