%% Error Calculations Second and Fourth
clc;clear;
%% Initializations
alpha = 1;

pol_signal_one = [1;-1;0]/sqrt(2);
pol_cum_signal_one = [1;1;0;-1;0;0]/2; %ground

% pol_signal_one = [1;1;1]/sqrt(3);
% pol_cum_signal_one = [1;1;1;1;1;1]/3;

signal_one_offset = 30*pi/180;
samples = 500;
Window = 100;    %size of window
Averaged_samples = 100;

signal_one_phase_est_second = zeros(1,Averaged_samples);
signal_one_mag_est_second = zeros(1,Averaged_samples);

signal_one_phase_est_fourth = zeros(1,Averaged_samples);
signal_one_mag_est_fourth = zeros(1,Averaged_samples);

SNR = zeros(1,Averaged_samples);
%% Matrix Calculations
% Implimenting a window. Esprit and SR techniques
Noise = linspace(0.1,10,Averaged_samples);
for Averaged_sample = 1:Averaged_samples;
    
    SNR(Averaged_sample) = 10*log10(1/(Noise(Averaged_sample))^2);
    for n = 1:samples
        %% Random Statistics Used for Both second and forth order Algorithms
        signal_one =  pol_signal_one*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        
        s1 = alpha*signal_one;
        s2 = alpha*exp(1i*signal_one_offset)*signal_one ;
        
        s1_Noise = s1 + Noise(Averaged_sample)*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        s2_Noise = s2 + Noise(Averaged_sample)*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        %% Second order Statistics
        p1_2 = [s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        p2_2 = [s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        R1_2 = p1_2*p1_2'/Window;
        R2_2 = p1_2*p2_2'/Window;
        
        A_2 = pinv(R1_2)*R2_2;
        
        [u_2,c_2] = eig(A_2);
        
        sort_pol_2 = abs(pol_signal_one'*u_2);
        [~,kk_2] = sort(sort_pol_2,'descend');
        signal_one_phase_est_second(Averaged_sample) = signal_one_phase_est_second(Averaged_sample) + ((angle(c_2(kk_2(1),kk_2(1)))+signal_one_offset)^2)/samples;
        signal_one_mag_est_second(Averaged_sample) = signal_one_mag_est_second(Averaged_sample) + (abs(c_2(kk_2(1),kk_2(1))))/samples;
        
        %% Fourth Order Statistics
        p1_4 = [s1_Noise(1,:).*s1_Noise(1,:)
            s1_Noise(2,:).*s1_Noise(2,:)
            s1_Noise(3,:).*s1_Noise(3,:)
            s1_Noise(1,:).*s1_Noise(2,:)
            s1_Noise(1,:).*s1_Noise(3,:)
            s1_Noise(2,:).*s1_Noise(3,:)];
        
        p2_4 = [s2_Noise(1,:).*s2_Noise(1,:)
            s2_Noise(2,:).*s2_Noise(2,:)
            s2_Noise(3,:).*s2_Noise(3,:)
            s2_Noise(1,:).*s2_Noise(2,:)
            s2_Noise(1,:).*s2_Noise(3,:)
            s2_Noise(2,:).*s2_Noise(3,:)];
        
        R1_4 = p1_4*p1_4'/Window;
        R2_4 = p1_4*p2_4'/Window;
        
        A_4 = pinv(R1_4)*R2_4;
        
        [u_4,c_4] = eig(A_4);
        
        sort_pol_4 = abs(pol_cum_signal_one'*u_4);
        [~,kk_4] = sort(sort_pol_4,'descend');
        signal_one_phase_est_fourth(Averaged_sample) = signal_one_phase_est_fourth(Averaged_sample) + ((0.5*angle(c_4(kk_4(1),kk_4(1)))+signal_one_offset)^2)/samples;
        signal_one_mag_est_fourth(Averaged_sample) = signal_one_mag_est_fourth(Averaged_sample) + (abs(c_4(kk_4(1),kk_4(1))))/samples;
        
    end
end
signal_one_phase_est_second_mnsqrd = sqrt(signal_one_phase_est_second)*180/pi;
signal_one_mag_est_second_mnsqrd = sqrt(signal_one_mag_est_second)*180/pi;
signal_one_phase_est_fourth_mnsqrd = sqrt(signal_one_phase_est_fourth)*180/pi;
signal_one_mag_est_fourth_mnsqrd = sqrt(signal_one_mag_est_fourth)*180/pi;
%% Plotting Results

figure(1);title('Mean Square Error (degrees)')
plot(SNR,signal_one_phase_est_second_mnsqrd,'b');
hold on;
plot(SNR,signal_one_phase_est_fourth_mnsqrd,'g');
hold off;
xlabel('SNR');ylabel('Mean Square Error (degrees)')
legend('2nd Order Error','4rth Order Error','Location','northeast');

figure(2);title('Mean Square Error (dB)')
plot(SNR,log10(signal_one_phase_est_second_mnsqrd),'b');
hold on;
plot(SNR,log10(signal_one_phase_est_fourth_mnsqrd),'g');
hold off;
xlabel('SNR');ylabel('Mean Square Error (dB)')
legend('2nd Order Error','4rth Order Error','Location','northeast');