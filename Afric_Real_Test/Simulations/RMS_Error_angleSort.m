%% Error Calculations Second and Fourth
clc;clear;
%% Initializations
alpha = 1;
beta = 1;

Pol_ground = [1;-1;0]/sqrt(2);
Pol_vegitation = [1;1;1]/sqrt(3);

G_O = 20;
ground_offset = G_O*pi/180; % ground interferomitry offset
V_O = 25;
vegitation_offset = V_O*pi/180;    % veg interferomitry offset


Averaged_samples = 1000;
Window = 81;    %size of window
SNR_samples = 30;

ground_phase_est_second = zeros(1,SNR_samples);
ground_mag_est_second = zeros(1,SNR_samples);

ground_phase_est_fourth = zeros(1,SNR_samples);
ground_mag_est_fourth = zeros(1,SNR_samples);

SNR = zeros(1,SNR_samples);
%% Generic Loop Calculations

for SNR_sample = fliplr(1:SNR_samples);
    
    SNR(SNR_sample)=SNR_sample-20;
    Noise = (10^(-SNR(SNR_sample)/20))/sqrt(3);
    
    for sample = 1:Averaged_samples
        %% Random Statistics Used for Both second and forth order Algorithms
        
        g =  Pol_ground*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        v =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        
        s1 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        %% Second order Statistics
        S1_2 = [s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        S2_2 = [s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        R1_2 = S1_2*S1_2'/Window;
        R2_2 = S1_2*S2_2'/Window;
        
        [eigenvect_2,eigenval_2] = eig(pinv(R1_2)*R2_2);
        
        [~,srt_2] = sort(abs(ground_offset - abs(angle(diag(eigenval_2)))),'ascend');
        
        ground_phase_est_second(SNR_sample) = ground_phase_est_second(SNR_sample)...
            + ((ground_offset - abs(angle(eigenval_2(srt_2(1),srt_2(1)))))^2)/Averaged_samples;
        
        ground_mag_est_second(SNR_sample) = ground_mag_est_second(SNR_sample)...
            + (abs(eigenval_2(srt_2(1),srt_2(1))))/Averaged_samples;
        
        %% Fourth Order Statistics
        S1_4 = [s1_Noise(1,:).*s1_Noise(1,:)
            s1_Noise(2,:).*s1_Noise(2,:)
            s1_Noise(3,:).*s1_Noise(3,:)
            s1_Noise(1,:).*s1_Noise(2,:)
            s1_Noise(1,:).*s1_Noise(3,:)
            s1_Noise(2,:).*s1_Noise(3,:)];
        
        S2_4 = [s2_Noise(1,:).*s2_Noise(1,:)
            s2_Noise(2,:).*s2_Noise(2,:)
            s2_Noise(3,:).*s2_Noise(3,:)
            s2_Noise(1,:).*s2_Noise(2,:)
            s2_Noise(1,:).*s2_Noise(3,:)
            s2_Noise(2,:).*s2_Noise(3,:)];
        
        R1_4 = S1_4*S1_4'/Window;
        R2_4 = S1_4*S2_4'/Window;
        
        [eigenvec_4,eigenval_4] = eig(pinv(R1_4)*R2_4);
               
        [~,srt_4]=sort(abs(ground_offset - abs(0.5*angle(diag(eigenval_4)))),'ascend');
        
        ground_phase_est_fourth(SNR_sample) = ground_phase_est_fourth(SNR_sample)...
            + ((ground_offset - abs(0.5*angle(eigenval_4(srt_4(1),srt_4(1)))))^2)/Averaged_samples;
        
        ground_mag_est_fourth(SNR_sample) = ground_mag_est_fourth(SNR_sample)...
            + (abs(eigenval_4(srt_4(1),srt_4(1))))/Averaged_samples;
        
    end
end
ground_phase_est_second_rmnsqrd = sqrt(ground_phase_est_second)*180/pi;
ground_phase_est_fourth_mnsqrd = sqrt(ground_phase_est_fourth)*180/pi;

%% Plotting Results
%load('RMS_Error10Ksamples.mat')
figure(1);title('RMS Error (dB) Vs SNR')
plot(SNR,10*log10(ground_phase_est_second_rmnsqrd),'b');
hold on;
plot(SNR,10*log10(ground_phase_est_fourth_mnsqrd),'black');
hold off;
xlabel('SNR');ylabel('RMS Error(dB)')
legend('2nd Order Error','4rth Order Error','Location','northeast');

figure(2);title('RMS Error vs Coherence')
plot(ground_mag_est_second,ground_phase_est_second_rmnsqrd,'b');
hold on;
plot(ground_mag_est_fourth,ground_phase_est_fourth_mnsqrd,'black');
hold off;
ylabel('RMS Error (degrees)');xlabel('Coherence')
legend('2nd Order Error Vs Coherence','4rth Order Error Vs Coherence','Location','northeast');
