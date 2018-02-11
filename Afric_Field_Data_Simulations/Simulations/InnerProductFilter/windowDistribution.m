%% Error Calculations Second and Fourth
clc;clear;
%% Initializations
alpha = 1;
eye_optimal = 0.3737;

pol_signal_one = [1;-1;0]./sqrt(2);
pol_cum_signal_one = [1;1;0;-1;0;0]./sqrt(3); %ground

signal_one_offset = 30*pi/180;

Averaged_samples = 1000;
window_dist = 100;
win = 10:window_dist;

phase_dist_second_10 = zeros(1,window_dist);
mag_dist_second_10 = zeros(1,window_dist);

phase_dist_second_0 = zeros(1,window_dist);
mag_dist_second_0 = zeros(1,window_dist);

phase_dist_second_n10 = zeros(1,window_dist);
mag_dist_second_n10 = zeros(1,window_dist);

phase_dist_fourth_10 = zeros(1,window_dist);
mag_dist_fourth_10 = zeros(1,window_dist);

phase_dist_fourth_0 = zeros(1,window_dist);
mag_dist_fourth_0 = zeros(1,Averaged_samples);

phase_dist_fourth_n10 = zeros(1,window_dist);
mag_dist_fourth_n10 = zeros(1,window_dist);

SNR = [-10 ,0, 10];

%% Generic Loop Calculations

for SNR_sample = 1:length(SNR);
    
    Noise = (10^(-SNR(SNR_sample)/20))/sqrt(3);
    for window = win
        for sample = 1:Averaged_samples
            %% Random Statistics Used for Both second and forth order Algorithms
            
            signal_one =  pol_signal_one*(sqrt(-2*log(1-rand(1,window))).*exp(1i*2*pi*rand(1,window)));
            
            s1 = alpha*signal_one;
            s2 = alpha*exp(1i*signal_one_offset)*signal_one ;
            
            s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(3,window))).*exp(1i*2*pi*rand(3,window));
            s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,window))).*exp(1i*2*pi*rand(3,window));
            %% Second order Statistics
            S1_2 = [s1_Noise(1,:)
                s1_Noise(2,:)
                s1_Noise(3,:)];
            
            S2_2 = [s2_Noise(1,:)
                s2_Noise(2,:)
                s2_Noise(3,:)];
            
            R1_2 = S1_2*S1_2'/window;
            R2_2 = S1_2*S2_2'/window;
            
            [eigenvect_2,eigenval_2] = eig(pinv(R1_2 + eye_optimal*eye(3))*R2_2);
            
            polarfilter_2 = abs(pol_signal_one'*eigenvect_2);
            [~,srt_2] = sort(polarfilter_2,'descend');
            if (SNR(SNR_sample) == 10)
                phase_dist_second_10(window) = phase_dist_second_10(window) + ((signal_one_offset - abs(angle(eigenval_2(srt_2(1),srt_2(1)))))^2)/Averaged_samples;
                mag_dist_second_10(window) =  mag_dist_second_10(window) + (abs(eigenval_2(srt_2(1),srt_2(1))))/Averaged_samples;
                
            elseif (SNR(SNR_sample) == 0)
                phase_dist_second_0(window) = phase_dist_second_0(window) + ((signal_one_offset - abs(angle(eigenval_2(srt_2(1),srt_2(1)))))^2)/Averaged_samples;
                mag_dist_second_0(window) = mag_dist_second_0(window) + (abs(eigenval_2(srt_2(1),srt_2(1))))/Averaged_samples;
                
                
            elseif (SNR(SNR_sample) == -10)
                phase_dist_second_n10(window) = phase_dist_second_n10(window) + ((signal_one_offset - abs(angle(eigenval_2(srt_2(1),srt_2(1)))))^2)/Averaged_samples;
                mag_dist_second_n10(window) = mag_dist_second_n10(window) + (abs(eigenval_2(srt_2(1),srt_2(1))))/Averaged_samples;
                
            end
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
            
            R1_4 = S1_4*S1_4'/window;
            R2_4 = S1_4*S2_4'/window;
            
            [eigenvec_4,eigenval_4] = eig(pinv(R1_4 + eye_optimal*eye(6))*R2_4);
            
            polarfilter_4 = abs(pol_cum_signal_one'*eigenvec_4);
            [~,srt_4] = sort(polarfilter_4,'descend');
            if (SNR(SNR_sample) == 10)
                phase_dist_fourth_10(window) = phase_dist_fourth_10(window) + ((signal_one_offset - abs(0.5*angle(eigenval_4(srt_2(1),srt_2(1)))))^2)/Averaged_samples;
                mag_dist_fourth_10(window) = mag_dist_fourth_10(window) + (abs(eigenval_4(srt_2(1),srt_2(1))))/Averaged_samples;
                
            elseif (SNR(SNR_sample) == 0)
                phase_dist_fourth_0(window) = phase_dist_fourth_0(window) + ((signal_one_offset - abs(0.5*angle(eigenval_4(srt_2(1),srt_2(1)))))^2)/Averaged_samples;
                mag_dist_fourth_0(window) =  mag_dist_fourth_0(window) + (abs(eigenval_4(srt_2(1),srt_2(1))))/Averaged_samples;
                
                
            elseif (SNR(SNR_sample) == -10)
                phase_dist_fourth_n10(window) =  phase_dist_fourth_n10(window) + ((signal_one_offset - abs(0.5*angle(eigenval_4(srt_2(1),srt_2(1)))))^2)/Averaged_samples;
                mag_dist_fourth_n10(window) =  mag_dist_fourth_n10(window) + (abs(eigenval_4(srt_2(1),srt_2(1))))/Averaged_samples;
                
            end
        end
    end
end
phase_est_second_10_rms = sqrt(phase_dist_second_10)*180/pi;
phase_est_second_0_rms = sqrt(phase_dist_second_0)*180/pi;
phase_est_second_n10_rms = sqrt(phase_dist_second_n10)*180/pi;

phase_est_fourth_10_rms = sqrt(phase_dist_fourth_10)*180/pi;
phase_est_fourth_0_rms = sqrt(phase_dist_fourth_0)*180/pi;
phase_est_fourth_n10_rms = sqrt(phase_dist_fourth_n10)*180/pi;

%% Plotting Results
 load('windowDistribution10Ksamples.mat');

figure(1);
subplot(2,1,1);
plot(win,phase_est_second_10_rms);
title('window dist second Ord rms error (degrees) SNR 10 (dB)');
xlabel('Window Size');ylabel('Averave RMS error (degrees)');

% figure(4);
subplot(2,1,2);
plot(win,phase_est_fourth_10_rms);
title('window dist fourth Ord rms error (degrees) SNR 10 (dB)');
xlabel('Window Size');ylabel('Averave RMS error (degrees)');

figure(2);
subplot(2,1,1);
plot(win,phase_est_second_0_rms);
title('window dist second Ord rms error (degrees) SNR 0 (dB)');
xlabel('Window Size');ylabel('Averave RMS error (degrees)');

% figure(5);
subplot(2,1,2);
plot(win,phase_est_fourth_0_rms);
title('window dist fourth Ord rms error (degrees) SNR 0 (dB)');
xlabel('Window Size');ylabel('Averave RMS error (degrees)');


figure(3);
subplot(2,1,1);
plot(win,phase_est_second_n10_rms);
title('window dist second Ord rms error (degrees) SNR -10 (dB)');
xlabel('Window Size');ylabel('Averave RMS error (degrees)');

% figure(6);
subplot(2,1,2);
plot(win,phase_est_fourth_n10_rms);
title('window dist fourth Ord rms error (degrees) SNR -10 (dB)');
xlabel('Window Size');ylabel('Averave RMS error (degrees)');