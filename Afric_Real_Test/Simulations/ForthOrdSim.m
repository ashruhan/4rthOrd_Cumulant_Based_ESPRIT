%% Afric_Sim_test_ESPRIT
clc;clear;
%% Initializations

Averaged_samples = 1;
Window = 81;    %size of window
SNR_samples = 30;

alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor
G_O = 20;
ground_offset = G_O*pi/180; % ground interferomitry offset
V_O = 60;
vegitation_offset = V_O*pi/180;    % veg interferomitry offset
G_O_array = ones(SNR_samples,1)*G_O;
V_O_array = ones(SNR_samples,1)*V_O;

Pol_ground = [1;-1;0]/sqrt(2);
Pol_ground_4 = [1;1;0;-1;0;0]/sqrt(3); %ground
Pol_vegitation = [1;1;1]/sqrt(3);
Pol_vegitation_4 = [1;1;1;1;1;1]/sqrt(6); %vegitation

ground_phase_4 = zeros(SNR_samples,1);
ground_mag_4 = zeros(SNR_samples,1);
ground_phase_rmse_4 = zeros(SNR_samples,1);

vegitation_phase_4 = zeros(SNR_samples,1);
vegitation_mag_4 = zeros(SNR_samples,1);
vegitation_phase_rmse_4 = zeros(SNR_samples,1);

SNR = zeros(1,SNR_samples);
%% Matrix Construction
for Averaged_sample = 1:SNR_samples;
    
    SNR(Averaged_sample)=Averaged_sample-20;
    Noise = (10^(-SNR(Averaged_sample)/20))/sqrt(3);
    
    for unusedvariable = 1:Averaged_samples
        
        g =  Pol_ground*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        v =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        
        polfilter_4 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;
        
        s1_Noise = polfilter_4 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
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
        
        polfilter_4 = abs(Pol_ground_4'*eigenvec_4);
        [~,srt_4] = sort(polfilter_4,'descend');
        ground_phase_4(Averaged_sample) = ground_phase_4(Averaged_sample) + 0.5*angle(eigenval_4(srt_4(1),srt_4(1)))/Averaged_samples;
        ground_mag_4(Averaged_sample) = ground_mag_4(Averaged_sample) + abs(eigenval_4(srt_4(1),srt_4(1)))/Averaged_samples;
        ground_phase_rmse_4(Averaged_sample) =  ground_phase_rmse_4(Averaged_sample) + ((ground_offset + 0.5*angle(eigenval_4(srt_4(1),srt_4(1))))^2)/Averaged_samples;
        
        polfilter_4 = abs(Pol_vegitation_4'*eigenvec_4);
        [~,srt_4] = sort(polfilter_4,'descend');
        vegitation_phase_4(Averaged_sample) = vegitation_phase_4(Averaged_sample) + 0.5*angle(eigenval_4(srt_4(1),srt_4(1)))/Averaged_samples;
        vegitation_mag_4(Averaged_sample) = vegitation_mag_4(Averaged_sample) + abs(eigenval_4(srt_4(1),srt_4(1)))/Averaged_samples;
        vegitation_phase_rmse_4(Averaged_sample) =  vegitation_phase_rmse_4(Averaged_sample) + ((vegitation_offset + 0.5*angle(eigenval_4(srt_4(1),srt_4(1))))^2)/Averaged_samples;
        
    end
end
ground_phase_rmse_4 = sqrt(ground_phase_rmse_4)*180/pi;
vegitation_phase_rmse_4 = sqrt(vegitation_phase_rmse_4)*180/pi;
%% Plotting Results

figure(1);
title('4rth Order ESPRIT Interferometric Phases');
xlabel('SNR dB');ylabel('Interferometric Phase (Degrees)');
hold on;
plot(SNR,G_O_array,'b+');
plot(SNR,-1*ground_phase_4*180/pi,'bo');
plot(SNR,V_O_array,'g+');
plot(SNR,-1*vegitation_phase_4*180/pi,'go');
legend('4rth Order Ground','4rth Order Ground Estimated','2rth Order Vegitation','4rth Order Vegitaion Estimated','Location','southeast')
axis([-20, 10, -25 100]);
hold off

figure(2)
title('4rth Order ESPRIT RMS Error')
xlabel('SNR (dB)');ylabel('RMS Error (Degrees)');
hold on;
plot(SNR,ground_phase_rmse_4,'bo');
plot(SNR,vegitation_phase_rmse_4,'go');
axis([-20, 10, 0, 100]);
hold off;
legend('4rth Order Ground RMS Error','4rth Order Vegitation RMS Error','Location','northeast');

figure(3)
title('4rth Order ESPRIT RMS Error vs Coherence')
xlabel('Coherence');ylabel('RMS Error (Degrees)');
hold on;
plot(ground_mag_4,ground_phase_rmse_4,'bo');
plot(vegitation_mag_4,vegitation_phase_rmse_4,'go');
axis([0, 1, 0, 100]);
hold off;
legend('4rth Order Ground RMS Error','4rth Order Vegitation RMS Error','Location','northeast');

figure(4)
title('4rth Order ESPRIT Coherance');
xlabel('SNR (dB)');ylabel('Maginitude')
hold on;
plot(SNR,ground_mag_4,'bo');
plot(SNR,vegitation_mag_4,'go');
legend('4rth Order Ground Coherancee','4rth Order Vegitation Coherancee','Location','southeast');
axis([-20, 10, 0 1]);
hold off;