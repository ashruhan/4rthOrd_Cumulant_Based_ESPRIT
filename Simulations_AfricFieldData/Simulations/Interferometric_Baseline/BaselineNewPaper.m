
%% ESPRIT_baseline
clear;clc;
%% Initializations
signal_one = 1;
signal_one_angle = 25*pi/180;
delx = 0.5; % Delta x of baseline
baseline = 2*delx; %Physical baseline of two sample arrays

signal_one_phase = sin(signal_one_angle)*2*pi*baseline; %Used in the sorting routine
array_size = 6; % the amount of ground samples

samples = 100;
window = samples; %Averaged_Signal_samples in matrix
SNR_samples = 30;

est_signal_one_rms = zeros(1,SNR_samples);
est_sig_one_angle = zeros(1,SNR_samples);
est_sig_one_coherence = zeros(SNR_samples,1);
SNR = zeros(1,SNR_samples);

%% Start of Algorithm
for SNR_sample = 1:SNR_samples;
    SNR(SNR_sample) = SNR_sample-20;
    Noise = (10^(-SNR(SNR_sample)/20))/sqrt(3);
    for unusedVariable = 1:samples;
        %% Setting up Random Phi that differ across samples
        phi_one = 1i*2*pi*rand(1,window);
        pos = [1:array_size]'; %#ok<NBRAK>
        distance = 2.*pi.*pos.*delx;  % element separation of 1/2 wavelength
        
        %% Element Array size with Signal
        
        X = signal_one.*(exp(-1i*distance*sin(signal_one_angle))*exp(phi_one))... %Ground
            + Noise.*sqrt(-2.*log(1-rand(array_size,window))).*exp(1i*2*pi*rand(array_size,window));  %Noise
        
        %% Setting up Baselines
        % Remember to change the Initial Variablev delx if Baseline changes
        
        S1 = [X(2,:);X(3,:)];
        S2 = [X(4,:);X(5,:)];
        
        %% ESPRIT Algorithm
        R1 = S1*S1'/window;
        R2 = S1*S2'/window;
        
        [u,uv] = eig(pinv(R1)*R2);
        
        [~,kk]=sort(abs(angle(diag(uv))-signal_one_phase),'ascend');
        est_sig_one_angle(SNR_sample) = est_sig_one_angle(SNR_sample)...
            + (asin(angle(uv(kk(1),kk(1)))/(2*pi*baseline)))/samples;
        
        est_signal_one_rms(SNR_sample) = est_signal_one_rms(SNR_sample)...
            + ((signal_one_angle - (asin(angle(uv(kk(1),kk(1)))/(2*pi*baseline))))^2)/samples;
        
        est_sig_one_coherence(SNR_sample) = est_sig_one_coherence(SNR_sample) + abs(uv(kk(1),kk(1)))/samples;
        
    end
end
est_signal_one_rms = sqrt(est_signal_one_rms)*180/pi;
%% Ploting Results

figure(1);
plot(SNR,est_sig_one_angle*180/pi,'bo');
hold on;
plot(SNR,ones(1,SNR_samples)*signal_one_angle*180/pi,'b+');
title('Direction of Arrival Resolution');
xlabel('SNR (dB)');ylabel('Angle from Referance Plane (Degrees)');
legend('Signal One Estimated','Signal One Actual','Location','east')
hold off;

figure(2)
title('DOA Coherance');
xlabel('SNR (dB)');ylabel('Maginitude')
hold on;
plot(SNR,est_sig_one_coherence,'ro');
legend('Signal Coherance','Location','east');
hold off;

figure(3)
plot(SNR,(est_signal_one_rms),'bo')
title('RMS Error vs SNR');
xlabel('SNR (dB)');ylabel('RMS Error (Degrees)');

figure(4)
plot(est_sig_one_coherence,(est_signal_one_rms),'bo')
title('RMS Error vs Coherence');
xlabel('Coherence');ylabel('RMS Error (Degrees)');
