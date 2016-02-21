%% Error Calculations Second and Fourth
clc;clear;
%% Initializations
alpha = 1;

pol_signal_one = [1;-1;0]./sqrt(2);
pol_cum_signal_one = [1;1;0;-1;0;0]./sqrt(3); %ground

S_O = 30;
signal_one_offset = S_O*pi/180;

Averaged_samples = 1000;
Window = 49;    %size of window


phase_dist_second_10 = zeros(1,Averaged_samples);
mag_dist_second_10 = zeros(1,Averaged_samples);

phase_dist_second_0 = zeros(1,Averaged_samples);
mag_dist_second_0 = zeros(1,Averaged_samples);

phase_dist_second_n10 = zeros(1,Averaged_samples);
mag_dist_second_n10 = zeros(1,Averaged_samples);

phase_dist_fourth_10 = zeros(1,Averaged_samples);
mag_dist_fourth_10 = zeros(1,Averaged_samples);

phase_dist_fourth_0 = zeros(1,Averaged_samples);
mag_dist_fourth_0 = zeros(1,Averaged_samples);

phase_dist_fourth_n10 = zeros(1,Averaged_samples);
mag_dist_fourth_n10 = zeros(1,Averaged_samples);

SNR = [-10 ,0, 10];
%% Generic Loop Calculations

for SNR_sample = 1:length(SNR);
    
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
        if (SNR(SNR_sample) == 10)
            phase_dist_second_10(sample) = ((signal_one_offset + angle(eigenval_2(srt_2(1),srt_2(1)))));
            mag_dist_second_10(sample) =  (abs(eigenval_2(srt_2(1),srt_2(1))));
            
        elseif (SNR(SNR_sample) == 0)
            phase_dist_second_0(sample) = ((signal_one_offset + angle(eigenval_2(srt_2(1),srt_2(1)))));
            mag_dist_second_0(sample) = (abs(eigenval_2(srt_2(1),srt_2(1))));
            
            
        elseif (SNR(SNR_sample) == -10)
            phase_dist_second_n10(sample) = ((signal_one_offset + angle(eigenval_2(srt_2(1),srt_2(1)))));
            mag_dist_second_n10(sample) = (abs(eigenval_2(srt_2(1),srt_2(1))));
            
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
        
        R1_4 = S1_4*S1_4'/Window;
        R2_4 = S1_4*S2_4'/Window;
        
        [eigenvec_4,eigenval_4] = eig(pinv(R1_4)*R2_4);
        
        polarfilter_4 = abs(pol_cum_signal_one'*eigenvec_4);
        [~,srt_4] = sort(polarfilter_4,'descend');
        if (SNR(SNR_sample) == 10)
            phase_dist_fourth_10(sample) = ((signal_one_offset + 0.5*angle(eigenval_4(srt_2(1),srt_2(1)))));
            mag_dist_fourth_10(sample) = (abs(eigenval_4(srt_2(1),srt_2(1))));
            
        elseif (SNR(SNR_sample) == 0)
            phase_dist_fourth_0(sample) = ((signal_one_offset + 0.5*angle(eigenval_4(srt_2(1),srt_2(1)))));
            mag_dist_fourth_0(sample) =  (abs(eigenval_4(srt_2(1),srt_2(1))));
            
            
        elseif (SNR(SNR_sample) == -10)
            phase_dist_fourth_n10(sample) = ((signal_one_offset + 0.5*angle(eigenval_4(srt_2(1),srt_2(1)))));
            mag_dist_fourth_n10(sample) = (abs(eigenval_4(srt_2(1),srt_2(1))));
            
        end
    end
end
phase_dist_second_10 = (phase_dist_second_10)*180/pi;
phase_dist_second_0 = (phase_dist_second_0)*180/pi;
phase_dist_second_n10 = (phase_dist_second_n10)*180/pi;

phase_dist_fourth_10 = (phase_dist_fourth_10)*180/pi;
phase_dist_fourth_0 = (phase_dist_fourth_0)*180/pi;
phase_dist_fourth_n10 = (phase_dist_fourth_n10)*180/pi;

%% Plotting Results
bins = 50;
xmin = S_O-180;
xmax = S_O+180;
ymin = 0;
ymax = 300;

figure(1);
subplot(2,1,1);
hist(phase_dist_second_10,bins);
title('Phase Error Distribution SNR 10 (dB)');
xlabel('Degrees of Error');
legend('Second Order');
axis([xmin,xmax,ymin,ymax]);

% figure(4);
subplot(2,1,2);
hist(phase_dist_fourth_10,bins);
title('Phase Error Distribution SNR 10 (dB)');
xlabel('Degrees of Error');
legend('Fourth Order');
axis([xmin,xmax,ymin,ymax]);

ymax = 200;
figure(2);
subplot(2,1,1);
hist(phase_dist_second_0,bins);
title('Phase Error Distribution SNR 0 (dB)');
xlabel('Degrees of Error');
legend('Second Order');
axis([xmin,xmax,ymin,ymax]);

% figure(5);
subplot(2,1,2);
hist(phase_dist_fourth_0,bins);
title('Phase Error Distribution SNR 0 (dB)');
xlabel('Degrees of Error');
legend('Fourth Order');
axis([xmin,xmax,ymin,ymax]);

ymax = 150;
figure(3);
subplot(2,1,1);
hist(phase_dist_second_n10,bins);
title('Phase Error Distribution SNR -10 (dB)');
xlabel('Degrees of Error');
legend('Second Order');
axis([xmin,xmax,ymin,ymax]);

% figure(6);
subplot(2,1,2);
hist(phase_dist_fourth_n10,bins);
title('Phase Error Distribution SNR -10 (dB)');
xlabel('Degrees of Error');
legend('Fourth Order');
axis([xmin,xmax,ymin,ymax]);