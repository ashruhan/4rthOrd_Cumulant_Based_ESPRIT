%% Afric_Sim_test_ESPRIT
clc;clear;
%% Initializations
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor
%
Pol_ground = [1;-1;0]/sqrt(2);
Pol_Cum_ground = [1;1;0;-1;0;0]/sqrt(3); %ground
Pol_vegitation = [1;1;1]/sqrt(3);
Pol_Cum_vegitation = [1;1;1;1;1;1]/sqrt(6); %vegitation
G_O = 30;
ground_offset = G_O*pi/180; % ground interferomitry offset
V_O = 50;
vegitation_offset = V_O*pi/180;    % veg interferomitry offset
%
samples = 100;
Window = 49;    %size of window
Averaged_samples = 30;
%
ground_phase_cum=zeros(Averaged_samples,1); vegitation_phase_cum=zeros(Averaged_samples,1);
ground_mag_cum=zeros(Averaged_samples,1); vegitation_mag_cum=zeros(Averaged_samples,1);
ground_phase_esp=zeros(Averaged_samples,1); vegitation_phase_esp=zeros(Averaged_samples,1);
ground_mag_esp=zeros(Averaged_samples,1); vegitation_mag_esp=zeros(Averaged_samples,1);
%
SNR = zeros(1,Averaged_samples);
%% Matrix Construction
for Averaged_sample = 1:Averaged_samples;
    
    SNR(Averaged_sample)=Averaged_sample-10;
    Noise = (10^(-SNR(Averaged_sample)/20))/sqrt(3);
    
    for unusedvariable = 1:samples
        
        g =  Pol_ground*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        v =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        
        s1 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
 
 
%% Cumulant Algorithm        
        p1=[s1_Noise(1,:).*s1_Noise(1,:)
            s1_Noise(2,:).*s1_Noise(2,:)
            s1_Noise(3,:).*s1_Noise(3,:)
            s1_Noise(1,:).*s1_Noise(2,:)
            s1_Noise(1,:).*s1_Noise(3,:)
            s1_Noise(2,:).*s1_Noise(3,:)];
        
        p2=[s2_Noise(1,:).*s2_Noise(1,:)
            s2_Noise(2,:).*s2_Noise(2,:)
            s2_Noise(3,:).*s2_Noise(3,:)
            s2_Noise(1,:).*s2_Noise(2,:)
            s2_Noise(1,:).*s2_Noise(3,:)
            s2_Noise(2,:).*s2_Noise(3,:)];
        
        R11 = p1*p1'/Window;
        R22 = p1*p2'/Window;
        
        A = pinv(R11)*R22;
        
        [u,c] = eig(A);
        
        s1 = abs(Pol_Cum_ground'*u);
        [~,kk] = sort(s1,'descend');
        ground_phase_cum(Averaged_sample) = ground_phase_cum(Averaged_sample) + 0.5*angle(c(kk(1),kk(1)))/samples;
        ground_mag_cum(Averaged_sample) = ground_mag_cum(Averaged_sample) + abs(c(kk(1),kk(1)))/samples;
        
        s2 = abs(Pol_Cum_vegitation'*u);
        [~,kk] = sort(s2,'descend');
        vegitation_phase_cum(Averaged_sample) = vegitation_phase_cum(Averaged_sample) + 0.5*angle(c(kk(1),kk(1)))/samples;
        vegitation_mag_cum(Averaged_sample) = vegitation_mag_cum(Averaged_sample) + abs(c(kk(1),kk(1)))/samples;
        %% ESPRIT Algorithm
        p1 = [s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        p2 = [s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        R1 = p1*p1'/Window;
        R2 = p1*p2'/Window;
        
        A = pinv(R1)*R2;
        
        [u,c] = eig(A);
        
        s1 = abs(Pol_ground'*u);
        [~,kk] = sort(s1,'descend');
        ground_phase_esp(Averaged_sample) = ground_phase_esp(Averaged_sample) + angle(c(kk(1),kk(1)))/samples;
        ground_mag_esp(Averaged_sample) = ground_mag_esp(Averaged_sample) + abs(c(kk(1),kk(1)))/samples;
        
        s2 = abs(Pol_vegitation'*u);
        [~,kk] = sort(s2,'descend');
        vegitation_phase_esp(Averaged_sample) = vegitation_phase_esp(Averaged_sample) + angle(c(kk(1),kk(1)))/samples;
        vegitation_mag_esp(Averaged_sample) = vegitation_mag_esp(Averaged_sample) + abs(c(kk(1),kk(1)))/samples;
        
    end
end
%% Plotting Results
 
figure(1);
title('2nd and 4rth Order Modified ESPRIT Interferometric Phases');
xlabel('SNR dB');ylabel('Int Phase (Degrees)');
hold on;
plot(SNR,ground_phase_cum*180/pi,'bx');
plot(SNR,ground_phase_esp*180/pi,'bo');
plot(SNR,vegitation_phase_cum*180/pi,'gx');
plot(SNR,vegitation_phase_esp*180/pi,'go');
axis([-10 20 -60 -10]);
legend('4rth Order Ground','2nd Order Ground','4rth Order Vegitation','2nd Order Vegitaion','Location','northeast')
hold off
 
figure(2)
title('2nd and 4rth Order Modified ESPRIT Coherance');
xlabel('SNR (dB)');ylabel('Maginitude')
hold on;
plot(SNR,ground_mag_cum,'bx');
plot(SNR,ground_mag_esp,'bo');
plot(SNR,vegitation_mag_cum,'gx');
plot(SNR,vegitation_mag_esp,'go');
legend('4rth Order Ground Coherance','2nd Order Ground Coherancee','4rth Order Vegitation Coherance','2nd Order Vegitation Coherancee','Location','southeast');
hold off;
