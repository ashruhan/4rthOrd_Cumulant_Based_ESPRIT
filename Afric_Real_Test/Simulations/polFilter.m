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
Averaged_samples = 10;
Window = 49;    %size of window
SNR_samples = 30;
%
ground_phase_4 = zeros(SNR_samples,1);
ground_mag_4 = zeros(SNR_samples,1);
ground_phase_2 = zeros(SNR_samples,1);
ground_mag_2 = zeros(SNR_samples,1);

vegitation_phase_4 = zeros(SNR_samples,1);
vegitation_mag_4 = zeros(SNR_samples,1);
vegitation_phase_2 = zeros(SNR_samples,1);
vegitation_mag_2 = zeros(SNR_samples,1);
%
SNR = zeros(1,SNR_samples);
%% Matrix Construction
for Averaged_sample = 1:SNR_samples;
    
    SNR(Averaged_sample)=Averaged_sample-20;
    Noise = (10^(-SNR(Averaged_sample)/20))/sqrt(3);
    
    for unusedvariable = 1:Averaged_samples
        
        g =  Pol_ground*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        v =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        
        polfilter_2 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;
        
        s1_Noise = polfilter_2 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        
        
        %% Fourth Order ESPRIT
        S1_4=[s1_Noise(1,:).*s1_Noise(1,:)
            s1_Noise(2,:).*s1_Noise(2,:)
            s1_Noise(3,:).*s1_Noise(3,:)
            s1_Noise(1,:).*s1_Noise(2,:)
            s1_Noise(1,:).*s1_Noise(3,:)
            s1_Noise(2,:).*s1_Noise(3,:)];
        
        S2_4=[s2_Noise(1,:).*s2_Noise(1,:)
            s2_Noise(2,:).*s2_Noise(2,:)
            s2_Noise(3,:).*s2_Noise(3,:)
            s2_Noise(1,:).*s2_Noise(2,:)
            s2_Noise(1,:).*s2_Noise(3,:)
            s2_Noise(2,:).*s2_Noise(3,:)];
        
        R1_4 = S1_4*S1_4'/Window;
        R2_4 = S1_4*S2_4'/Window;
        
        [eigenvec_4,eigenval_4] = eig(pinv(R1_4)*R2_4);
        
        polfilter_4 = abs(Pol_Cum_ground'*eigenvec_4);
        [~,srt_4] = sort(polfilter_4,'descend');
        ground_phase_4(Averaged_sample) = ground_phase_4(Averaged_sample) + 0.5*angle(eigenval_4(srt_4(1),srt_4(1)))/Averaged_samples;
        ground_mag_4(Averaged_sample) = ground_mag_4(Averaged_sample) + abs(eigenval_4(srt_4(1),srt_4(1)))/Averaged_samples;
        
        polfilter_4 = abs(Pol_Cum_vegitation'*eigenvec_4);
        [~,srt_4] = sort(polfilter_4,'descend');
        vegitation_phase_4(Averaged_sample) = vegitation_phase_4(Averaged_sample) + 0.5*angle(eigenval_4(srt_4(1),srt_4(1)))/Averaged_samples;
        vegitation_mag_4(Averaged_sample) = vegitation_mag_4(Averaged_sample) + abs(eigenval_4(srt_4(1),srt_4(1)))/Averaged_samples;
        %% Second Order ESPRIT
        S1_2 = [s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        S2_2 = [s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        R1_2 = S1_2*S1_2'/Window;
        R2_2 = S1_2*S2_2'/Window;
        
        [eigenvec_2,eigenval_2] = eig(pinv(R1_2)*R2_2);
        
        polfilter_2 = abs(Pol_ground'*eigenvec_2);
        [~,srt_2] = sort(polfilter_2,'descend');
        ground_phase_2(Averaged_sample) = ground_phase_2(Averaged_sample) + angle(eigenval_2(srt_2(1),srt_2(1)))/Averaged_samples;
        ground_mag_2(Averaged_sample) = ground_mag_2(Averaged_sample) + abs(eigenval_2(srt_2(1),srt_2(1)))/Averaged_samples;
        
        polfilter_2 = abs(Pol_vegitation'*eigenvec_2);
        [~,srt_2] = sort(polfilter_2,'descend');
        vegitation_phase_2(Averaged_sample) = vegitation_phase_2(Averaged_sample) + angle(eigenval_2(srt_2(1),srt_2(1)))/Averaged_samples;
        vegitation_mag_2(Averaged_sample) = vegitation_mag_2(Averaged_sample) + abs(eigenval_2(srt_2(1),srt_2(1)))/Averaged_samples;
        
    end
end
%% Plotting Results

figure(1);
title('2nd and 4rth Order Modified ESPRIT Interferometric Phases');
xlabel('SNR dB');ylabel('Int Phase (Degrees)');
hold on;
plot(SNR,ground_phase_4*180/pi,'bx');
plot(SNR,ground_phase_2*180/pi,'bo');
plot(SNR,vegitation_phase_4*180/pi,'gx');
plot(SNR,vegitation_phase_2*180/pi,'go');
legend('4rth Order Ground','2nd Order Ground','4rth Order Vegitation','2nd Order Vegitaion','Location','northeast')
hold off

figure(2)
title('2nd and 4rth Order Modified ESPRIT Coherance');
xlabel('SNR (dB)');ylabel('Maginitude')
hold on;
plot(SNR,ground_mag_4,'bx');
plot(SNR,ground_mag_2,'bo');
plot(SNR,vegitation_mag_4,'gx');
plot(SNR,vegitation_mag_2,'go');
legend('4rth Order Ground Coherance','2nd Order Ground Coherancee','4rth Order Vegitation Coherance','2nd Order Vegitation Coherancee','Location','southeast');
hold off;