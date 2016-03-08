%% Afric_Sim_test_ESPRIT
clc;clear;
%% Initializations
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor
eye_optimal = 0.3737;

Pol_ground = [1;-1;0]./sqrt(2);
Pol_vegitation = [1;1;1]./sqrt(3);
G_O = 30;
ground_offset = G_O*pi/180; % ground interferomitry offset
V_O = 50;
vegitation_offset = V_O*pi/180;    % veg interferomitry offset

Averaged_samples = 10;
window_optimal = 81;    %size of window
SNR_samples = 30;
second_order = 3;

ground_phase_4=zeros(SNR_samples,1);
ground_mag_4=zeros(SNR_samples,1);
ground_phase_2=zeros(SNR_samples,1);
ground_mag_2=zeros(SNR_samples,1);

vegitation_phase_4=zeros(SNR_samples,1);
vegitation_mag_4=zeros(SNR_samples,1);
vegitation_phase_2=zeros(SNR_samples,1);
vegitation_mag_2=zeros(SNR_samples,1);

sort_ground_4 = zeros(2,1);
sort_vegetation_4 = zeros(2,1);

SNR = zeros(1,SNR_samples);
%% Matrix Construction
for SNR_sample = (1:SNR_samples);
    
    SNR(SNR_sample)=SNR_sample-20;
    Noise = (10^(-SNR(SNR_sample)/20))/sqrt(3);
    
    for unusedvariable = 1:Averaged_samples
        
        g =  Pol_ground*(sqrt(-2*log(1-rand(1,window_optimal))).*exp(1i*2*pi*rand(1,window_optimal)));
        v =  Pol_vegitation*(sqrt(-2*log(1-rand(1,window_optimal))).*exp(1i*2*pi*rand(1,window_optimal)));
        
        s1 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(second_order,window_optimal))).*exp(1i*2*pi*rand(second_order,window_optimal));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(second_order,window_optimal))).*exp(1i*2*pi*rand(second_order,window_optimal));
        
        
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
        
        R1_4 = S1_4*S1_4'/window_optimal;
        R2_4 = S1_4*S2_4'/window_optimal;
        
        [eigenvec_4,eigenval_4] = eig(pinv(R1_4 + eye_optimal*eye(6))*R2_4);
        
        [~,srtABS_4] = sort(abs(diag(eigenval_4)),'descend');
        
        Largest_eigenvals_4 = [eigenval_4(srtABS_4(1),srtABS_4(1))
                               eigenval_4(srtABS_4(2),srtABS_4(2))
                               eigenval_4(srtABS_4(3),srtABS_4(3))];
        
        [~,srtANGLE_4] = sort(abs(angle(Largest_eigenvals_4)),'descend');
        
        sort_ground_4(1) = (abs(eigenval_4(srtABS_4(srtANGLE_4(1)),srtABS_4(srtANGLE_4(1)))))^2....
                          *(abs(eigenvec_4(1,srtABS_4(srtANGLE_4(1))))^2....
                          + abs(eigenvec_4(2,srtABS_4(srtANGLE_4(1))))^2....
                          + abs(eigenvec_4(4,srtABS_4(srtANGLE_4(1))))^2);
        
        sort_ground_4(2) = abs(eigenval_4(srtABS_4(srtANGLE_4(3)),srtABS_4(srtANGLE_4(3))))^2....
                         *(abs(eigenvec_4(1,srtABS_4(srtANGLE_4(3))))^2....
                         + abs(eigenvec_4(2,srtABS_4(srtANGLE_4(3))))^2....
                         + abs(eigenvec_4(4,srtABS_4(srtANGLE_4(3))))^2);

        if sort_ground_4(1) >= sort_ground_4(2)
            ground_phase_4(SNR_sample) = ground_phase_4(SNR_sample)....
                + 0.5*angle(eigenval_4(srtABS_4(srtANGLE_4(1)),srtABS_4(srtANGLE_4(1))))/Averaged_samples;
        else
            ground_phase_4(SNR_sample) = ground_phase_4(SNR_sample)....
                + 0.5*angle(eigenval_4(srtABS_4(srtANGLE_4(3)),srtABS_4(srtANGLE_4(3))))/Averaged_samples;
        end
        
        sort_vegetation_4(1) = (abs(eigenval_4(srtABS_4(srtANGLE_4(1)),srtABS_4(srtANGLE_4(1)))))^2....
                              *(abs(eigenvec_4(3,srtABS_4(srtANGLE_4(1))))^2....
                              + abs(eigenvec_4(5,srtABS_4(srtANGLE_4(1))))^2....
                              + abs(eigenvec_4(6,srtABS_4(srtANGLE_4(1))))^2);
        
        sort_vegetation_4(2) = (abs(eigenval_4(srtABS_4(srtANGLE_4(3)),srtABS_4(srtANGLE_4(3)))))^2....
                              *(abs(eigenvec_4(3,srtABS_4(srtANGLE_4(3))))^2....
                              + abs(eigenvec_4(5,srtABS_4(srtANGLE_4(3))))^2....
                              + abs(eigenvec_4(6,srtABS_4(srtANGLE_4(3))))^2);

        if sort_vegetation_4(1) >= sort_vegetation_4(2)
            vegitation_phase_4(SNR_sample) = vegitation_phase_4(SNR_sample)....
                + 0.5*angle(eigenval_4(srtABS_4(srtANGLE_4(1)),srtABS_4(srtANGLE_4(1))))/Averaged_samples;
        else
            vegitation_phase_4(SNR_sample) = vegitation_phase_4(SNR_sample)....
                + 0.5*angle(eigenval_4(srtABS_4(srtANGLE_4(3)),srtABS_4(srtANGLE_4(3))))/Averaged_samples;
        end     
        %% Second Order ESPRIT
        S1_2 = [s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        S2_2 = [s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        R1_2 = S1_2*S1_2'/window_optimal;
        R2_2 = S1_2*S2_2'/window_optimal;
        
        [eigenvec_2,eigenval_2] = eig(pinv(R1_2 + eye_optimal*eye(3))*R2_2);
        
        [~,srt_2]=sort(abs(diag(eigenval_2)),'descend');
        
        Leig_copol = abs(eigenvec_2(1,srt_2(1)))^2 + abs(eigenvec_2(2,srt_2(1)))^2;
        SLeig_copol = abs(eigenvec_2(1,srt_2(2)))^2 + abs(eigenvec_2(2,srt_2(2)))^2;
        
        if (Leig_copol >= SLeig_copol)
            
            ground_phase_2(SNR_sample) = ground_phase_2(SNR_sample)....
                + angle(eigenval_2(srt_2(1),srt_2(1)))/Averaged_samples;
            
            ground_mag_2(SNR_sample) = ground_mag_2(SNR_sample)....
                + abs(eigenval_2(srt_2(1),srt_2(1)))/Averaged_samples;
            
            vegitation_phase_2(SNR_sample) = vegitation_phase_2(SNR_sample)....
                + angle(eigenval_2(srt_2(2),srt_2(2)))/Averaged_samples;
            
            vegitation_mag_2(SNR_sample) = vegitation_mag_2(SNR_sample)....
                + abs(eigenval_2(srt_2(2),srt_2(2)))/Averaged_samples;
            
        else
            
            ground_phase_2(SNR_sample) = ground_phase_2(SNR_sample)....
                + angle(eigenval_2(srt_2(2),srt_2(2)))/Averaged_samples;
            
            ground_mag_2(SNR_sample) = ground_mag_2(SNR_sample)....
                + abs(eigenval_2(srt_2(2),srt_2(2)))/Averaged_samples;
            
            vegitation_phase_2(SNR_sample) = vegitation_phase_2(SNR_sample)....
                + angle(eigenval_2(srt_2(1),srt_2(1)))/Averaged_samples;
            
            vegitation_mag_2(SNR_sample) = vegitation_mag_2(SNR_sample)....
                + abs(eigenval_2(srt_2(1),srt_2(1)))/Averaged_samples;
            
        end
        
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
plot(SNR,-V_O*ones(1,SNR_samples),'g');
plot(SNR,-G_O*ones(1,SNR_samples),'b');
legend('4rth Order Ground','2nd Order Ground','4rth Order Vegitation','2nd Order Vegitaion','Location','northeast')
hold off
% 
% figure(2);title('RMS Error (dB)')
% plot(SNR,(ground_phase_est_fourth_mnsqrd),'black');
% xlabel('SNR');ylabel('RMS Error(dB)')

% figure(2)
% title('2nd and 4rth Order Modified ESPRIT Coherance');
% xlabel('SNR (dB)');ylabel('Maginitude')
% hold on;
% plot(SNR,ground_mag_4,'bx');
% plot(SNR,ground_mag_2,'bo');
% plot(SNR,vegitation_mag_4,'gx');
% plot(SNR,vegitation_mag_2,'go');
% legend('4rth Order Ground Coherance','2nd Order Ground Coherancee','4rth Order Vegitation Coherance','2nd Order Vegitation Coherancee','Location','southeast');
% hold off;