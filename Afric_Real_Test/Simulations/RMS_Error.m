%% Error Calculations Second and Fourth
clc;clear;
%% Initializations
alpha = 1;
 
pol_signal_one = [1;-1;0]./sqrt(2);
pol_cum_signal_one = [1;1;0;-1;0;0]./sqrt(3); %ground

signal_one_offset = 30*pi/180;

Averaged_samples = 100;
Window = 49;    %size of window
SNR_samples = 30;
 
signal_one_phase_est_second = zeros(1,SNR_samples);
signal_one_mag_est_second = zeros(1,SNR_samples);
 
signal_one_phase_est_fourth = zeros(1,SNR_samples);
signal_one_mag_est_fourth = zeros(1,SNR_samples);
 
SNR = zeros(1,SNR_samples);
%% Generic Loop Calculations

for SNR_sample = 1:SNR_samples;
    
    SNR(SNR_sample)=SNR_sample-20;
    Noise = (10^(-SNR(SNR_sample)/20))/sqrt(3);
    
    for sample = 1:Averaged_samples
        %% Random Statistics Used for Both second and forth order Algorithms
        signal_one =  pol_signal_one*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        
        s1 = alpha*signal_one;
        s2 = alpha*exp(1i*signal_one_offset)*signal_one ;
        
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
        
        polarfilter_2 = abs(pol_signal_one'*eigenvect_2);
        [~,srt_2] = sort(polarfilter_2,'descend');
        signal_one_phase_est_second(SNR_sample) = signal_one_phase_est_second(SNR_sample) + ((signal_one_offset + angle(eigenval_2(srt_2(1),srt_2(1))))^2)/Averaged_samples;
        signal_one_mag_est_second(SNR_sample) = signal_one_mag_est_second(SNR_sample) + (abs(eigenval_2(srt_2(1),srt_2(1))))/Averaged_samples;
        
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
        
        polarfilter_4 = abs(pol_cum_signal_one'*eigenvec_4);
        [~,srt_4] = sort(polarfilter_4,'descend');
        signal_one_phase_est_fourth(SNR_sample) = signal_one_phase_est_fourth(SNR_sample) + ((signal_one_offset + 0.5*angle(eigenval_4(srt_4(1),srt_4(1))))^2)/Averaged_samples;
        signal_one_mag_est_fourth(SNR_sample) = signal_one_mag_est_fourth(SNR_sample) + (abs(eigenval_4(srt_4(1),srt_4(1))))/Averaged_samples;
        
    end
end
signal_one_phase_est_second_rmnsqrd = sqrt(signal_one_phase_est_second)*180/pi;
signal_one_phase_est_fourth_mnsqrd = sqrt(signal_one_phase_est_fourth)*180/pi;

%% Plotting Results
 
% figure(1);
% plot(SNR,log10(signal_one_phase_est_fourth_mnsqrd),'g');
% title(' Co and Cross Polarizations (dB)');
% xlabel('SNR (dB)');ylabel('Mean Square Error (dB)')
% legend('4rth Order Error','Location','northeast');
 
figure(1);
plot(SNR,signal_one_phase_est_second_rmnsqrd,'b');
hold on;
plot(SNR,signal_one_phase_est_fourth_mnsqrd,'black');
hold off;
xlabel('SNR');ylabel('RMS Error (degrees)')
legend('2nd Order Error','4rth Order Error','Location','northeast');
 
figure(2);title('RMS Error (dB)')
plot(SNR,10*log10(signal_one_phase_est_second_rmnsqrd),'b');
hold on;
plot(SNR,10*log10(signal_one_phase_est_fourth_mnsqrd),'black');
hold off;
xlabel('SNR');ylabel('RMS Error(dB)')
legend('2nd Order Error','4rth Order Error','Location','northeast');

figure(3);title('RMS Error vs Coherence')
plot(signal_one_mag_est_second,signal_one_phase_est_second_rmnsqrd,'b');
hold on;
plot(signal_one_mag_est_fourth,signal_one_phase_est_fourth_mnsqrd,'black');
hold off;
ylabel('RMS Error (degrees)');xlabel('Coherence')
legend('2nd Order Error Vs Coherence','4rth Order Error Vs Coherence','Location','northeast');
