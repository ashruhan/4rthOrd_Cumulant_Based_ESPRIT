%% Afric_Sim_test_ESPRIT
clc;clear;
%% Initializations

Averaged_samples = 10;
Window_optimal = 81;    %size of window
SNR_samples = 30;
eye_optimal = -0.022;

alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor
G_O = 20;
ground_offset = G_O*pi/180; % ground interferomitry offset
V_O = 60;
vegitation_offset = V_O*pi/180;    % veg interferomitry offset
G_O_array = ones(SNR_samples,1)*G_O;
V_O_array = ones(SNR_samples,1)*V_O;
Pol_ground = [1;1;0]/sqrt(2);
Pol_vegitation = [1;1;1]/sqrt(3);

ground_angle_2 = zeros(SNR_samples,1);
ground_abs_2 = zeros(SNR_samples,1);
ground_angle_rmse_2 = zeros(SNR_samples,1);

vegitation_angle_2 = zeros(SNR_samples,1);
vegitation_abs_2 = zeros(SNR_samples,1);
vegitation_angle_rmse_2 = zeros(SNR_samples,1);

SNR = zeros(1,SNR_samples);
%% Matrix Construction
for Averaged_sample = 1:SNR_samples;
    
    SNR(Averaged_sample)=Averaged_sample-15;
    Noise = (10^(-SNR(Averaged_sample)/20))/sqrt(3);
    
    for unusedvariable = 1:Averaged_samples
        
        g =  Pol_ground*exp(1i*2*pi*rand(1,Window_optimal));
        v =  Pol_vegitation*exp(1i*2*pi*rand(1,Window_optimal));
        
        s1 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(3,Window_optimal))).*exp(1i*2*pi*rand(3,Window_optimal));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,Window_optimal))).*exp(1i*2*pi*rand(3,Window_optimal));
       
        %% Second Order ESPRIT
        S1_2 = [s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        S2_2 = [s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        R1_2 = S1_2*S1_2'/Window_optimal;
        R2_2 = S1_2*S2_2'/Window_optimal;
        
        [eigenvec_2,eigenval_2] = eig(pinv(R1_2 + eye_optimal*eye(3))*R2_2,'nobalance');
        
        polfilter_2 = abs(Pol_ground'*eigenvec_2);
        [~,srt_2] = sort(polfilter_2,'descend');
        
        ground_angle_2(Averaged_sample) = ground_angle_2(Averaged_sample)...
            + angle(eigenval_2(srt_2(1),srt_2(1)))/Averaged_samples;
        
        ground_abs_2(Averaged_sample) = ground_abs_2(Averaged_sample)...
            + abs(eigenval_2(srt_2(1),srt_2(1)))/Averaged_samples;
        
        ground_angle_rmse_2(Averaged_sample) =  ground_angle_rmse_2(Averaged_sample)...
            + ((ground_offset + angle(eigenval_2(srt_2(1),srt_2(1))))^2)/Averaged_samples;
        
        polfilter_2 = abs(Pol_vegitation'*eigenvec_2);
        [~,srt_2] = sort(polfilter_2,'descend');
        
        vegitation_angle_2(Averaged_sample) = vegitation_angle_2(Averaged_sample)...
            + angle(eigenval_2(srt_2(1),srt_2(1)))/Averaged_samples;
        
        vegitation_abs_2(Averaged_sample) = vegitation_abs_2(Averaged_sample)...
            + abs(eigenval_2(srt_2(1),srt_2(1)))/Averaged_samples;
        
        vegitation_angle_rmse_2(Averaged_sample) =  vegitation_angle_rmse_2(Averaged_sample)...
            + ((vegitation_offset + angle(eigenval_2(srt_2(1),srt_2(1))))^2)/Averaged_samples;
        
    end
end
ground_angle_rmse_2 = sqrt(ground_angle_rmse_2)*180/pi;
vegitation_angle_rmse_2 = sqrt(vegitation_angle_rmse_2)*180/pi;
%% Plotting Results

figure(1);
title('2nd Order ESPRIT Interferometric Phases');
xlabel('SNR dB');ylabel('Interferometric Phase (Degrees)');
hold on;
plot(SNR,G_O_array,'b+');
plot(SNR,-1*ground_angle_2*180/pi,'bo');
plot(SNR,V_O_array,'g+');
plot(SNR,-1*vegitation_angle_2*180/pi,'go');
legend('2nd Order Ground','2nd Order Ground Estimated','2rth Order Vegitation','2nd Order Vegitaion Estimated','Location','southeast')
axis([-15, 15, -25 100]);
hold off

figure(2)
title('2nd Order ESPRIT RMS Error')
xlabel('SNR (dB)');ylabel('RMS Error (Degrees)');
hold on;
plot(SNR,ground_angle_rmse_2,'bo');
plot(SNR,vegitation_angle_rmse_2,'go');
axis([-15, 15, 0, 100]);
hold off;
legend('2nd Order Ground RMS Error','2nd Order Vegitation RMS Error','Location','northeast');

figure(3)
title('2nd Order ESPRIT RMS Error vs Coherence')
xlabel('Coherence');ylabel('RMS Error (Degrees)');
hold on;
plot(ground_abs_2,ground_angle_rmse_2,'bo');
plot(vegitation_abs_2,vegitation_angle_rmse_2,'go');
axis([0, 1, 0 100]);
hold off;
legend('2nd Order Ground','2nd Order Vegitation','Location','northeast');

figure(4)
title('2nd Order ESPRIT Coherance');
xlabel('SNR (dB)');ylabel('Maginitude')
hold on;
plot(SNR,ground_abs_2,'bo');
plot(SNR,vegitation_abs_2,'go');
legend('2nd Order Ground Coherancee','2nd Order Vegitation Coherancee','Location','southeast');
axis([-15, 15 , 0 1]);
hold off;