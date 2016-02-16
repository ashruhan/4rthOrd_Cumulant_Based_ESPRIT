%% Afric_Sim_test_ESPRIT
clc;clear;
%% Initializations
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor
%
Pol_ground = [1;1;0]/sqrt(2);
Pol_vegitation = [1;1;1]/sqrt(3);

G_O = 30;
ground_offset = G_O*pi/180; % ground interferomitry offset
V_O = 65;
vegitation_offset = V_O*pi/180;    % veg interferomitry offset
%
samplesperSNR = 200;
Window = 200;    %size of window
SNR_samples = 30;
%

ground_phase=zeros(SNR_samples,1); vegitation_phase=zeros(SNR_samples,1);
ground_mag=zeros(SNR_samples,1); vegitation_mag=zeros(SNR_samples,1);
%
SNR = zeros(1,SNR_samples);
%% Matrix Construction
for sample = 1:SNR_samples;
    
    SNR(sample)=sample-10;
    Noise = (10^(-SNR(sample)/20))/sqrt(3);
    
    for unusedvariable = 1:samplesperSNR
        
        ground_polNoise =  Pol_ground*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        vege_polNoise =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        
        s1 = alpha*ground_polNoise + beta*vege_polNoise;
        s2 = alpha*exp(1i*ground_offset)*ground_polNoise + beta*exp(1i*vegitation_offset)*vege_polNoise;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        
        %% ESPRIT Algorithm
        S1 = [s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        S2 = [s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        R1 = S1*S1'/Window;
        R2 = S1*S2'/Window;
        
        [eigenvector,eigenvalue] = eig(pinv(R1)*R2);
        
        [~,kk] = sort(angle(diag(eigenvalue)),'descend');
        %
        %         if (abs(eigenvector(1,kk(1)))^2+abs(eigenvector(2,kk(1)))^2)<...
        %                 (abs(eigenvector(1,kk(2)))^2+abs(eigenvector(2,kk(2)))^2)
        %
        %             ground_phase_esp(sample) = ground_phase_esp(sample) + angle(eigenvalue(kk(1),kk(1)))/samplesperSNR;
        %             ground_mag_esp(sample) = ground_mag_esp(sample) + abs(eigenvalue(kk(1),kk(1)))/samplesperSNR;
        %
        %             vegitation_phase_esp(sample) = vegitation_phase_esp(sample) + angle(eigenvalue(kk(2),kk(2)))/samplesperSNR;
        %             vegitation_mag_esp(sample) = vegitation_mag_esp(sample) + abs(eigenvalue(kk(2),kk(2)))/samplesperSNR;
        %
        %
        %         elseif (abs(eigenvector(1,kk(2)))^2+abs(eigenvector(2,kk(2)))^2)<...
        %                 (abs(eigenvector(1,kk(1)))^2+abs(eigenvector(2,kk(1)))^2)
        %
        %             ground_phase_esp(sample) = ground_phase_esp(sample) + angle(eigenvalue(kk(2),kk(2)))/samplesperSNR;
        %             ground_mag_esp(sample) = ground_mag_esp(sample) + abs(eigenvalue(kk(2),kk(2)))/samplesperSNR;
        %
        %             vegitation_phase_esp(sample) = vegitation_phase_esp(sample) + angle(eigenvalue(kk(1),kk(1)))/samplesperSNR;
        %             vegitation_mag_esp(sample) = vegitation_mag_esp(sample) + abs(eigenvalue(kk(1),kk(1)))/samplesperSNR;
        %         end
        %
        
        ground_phase(sample) = ground_phase(sample) + angle(eigenvalue(kk(1),kk(1)))/samplesperSNR;
        ground_mag(sample) = ground_mag(sample) + abs(eigenvalue(kk(1),kk(1)))/samplesperSNR;
        
        vegitation_phase(sample) = vegitation_phase(sample) + angle(eigenvalue(kk(2),kk(2)))/samplesperSNR;
        vegitation_mag(sample) = vegitation_mag(sample) + abs(eigenvalue(kk(2),kk(2)))/samplesperSNR;
        
    end
end
%% Plotting Results

figure(1);
title('2nd ESPRIT Interferometric Phases');
xlabel('SNR dB');ylabel('Int Phase (Degrees)');
hold on;
plot(SNR,ground_phase*180/pi,'b');
plot(SNR,vegitation_phase*180/pi,'g');
hold off

figure(2)
title('2nd Order ESPRIT Coherance');
xlabel('SNR (dB)');ylabel('Maginitude')
hold on;
plot(SNR,ground_mag,'b');
plot(SNR,vegitation_mag,'g');
hold off;
