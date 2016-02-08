%% Espirit_baseline
clear;clc;
%% Initializations
signal_one = 1;
% signal_two = 1;
 
signal_one_angle = 10*pi/180; 
% signal_two_angle = 2*pi/180; 
 
delx = 0.5; % Delta x of baseline
baseline = 4*delx; %Physical baseline of two sample arrays
 
signal_one_phase = sin(signal_one_angle)*2*pi*baseline; %Used in the sorting routine
% signal_two_phase = sin(signal_two_angle)*2*pi*baseline; %Used in the sorting routine
 
array_size = 6; % the amount of ground samples
samples = 100;
Averaged_Signal_samples = 20; %Averaged_Signal_samples in matrix
Averaged_samples = 30;
 
est_sig_one_angle = zeros(1,Averaged_samples);
est_sig_one_mag=zeros(Averaged_samples,1); 
% est_sig_two_angle = zeros(1,Noise_samples);
SNR = zeros(1,Averaged_samples);
%% Start of Algorithm
for Averaged_sample = 1:Averaged_samples;
    SNR(Averaged_sample) = Averaged_sample-10;
    Noise = (10^(-SNR(Averaged_sample)/20))/sqrt(3);
    for unusedVariable = 1:samples;
        %% Setting up Random Phi that differ across samples
        phi_one = 1i*2*pi*rand(1,Averaged_Signal_samples);
%         phi_two = 1i*2*pi*rand(1,Averaged_Signal_samples);
        pos = [1:array_size]'; %#ok<NBRAK>
        distance = 2.*pi.*pos.*delx;  % element separation of 1/2 wavelength
        
        %% Element Array size with Signal
        
        X = signal_one.*(exp(-1i*distance*sin(signal_one_angle))*exp(phi_one))... %Ground
            + Noise.*sqrt(-2.*log(1-rand(array_size,Averaged_Signal_samples))).*exp(1i*2*pi*rand(array_size,Averaged_Signal_samples));  %Noise
        %             + signal_two.*(exp(-1i*distance*sin(signal_two_angle))*exp(phi_two))...  %Vegetaion
 
        %% Setting up Baselines
        % Remember to change the Initial Variable delx if Baseline changes
        
%                 S1 = [X(1,:); X(2,:)]; S2 = [X(5,:);X(6,:)]; %Baseline = 4delx
        %         S1 = [X(2,:);X(3,:)]; S2 = [X(4,:);X(5,:)]; %Baseline = 2delx
        S1 = [X(3,:);X(4,:)]; S2 = [X(4,:);X(5,:)]; %Baseline = delx
        %         S1 = [X(2,:);X(3,:);X(4,:)]; S2 = [X(3,:);X(4,:);X(5,:)]; %Baseline = delx
        %         S1 = [X(1,:);X(2,:);X(3,:);X(4,:)]; S2 = [X(2,:);X(3,:);X(4,:);X(5,:)]; %Baseline = delx
        
        %% Espirit Algorithm
        R1 = S1*S1'/Averaged_samples;
        R2 = S1*S2'/Averaged_samples;
        A =  pinv(R1)*R2;
        
        %% Averaging over a matix then averaging that signal in a for loop proportional to the matix rows
        % The Smaller Angle ends up haveing The largest angle Eigen Value
        
        [u,uv] = eig(A);
        
        [~,kk]=sort(abs(angle(diag(uv))-signal_one_phase),'ascend');
        est_sig_one_angle(Averaged_sample) = est_sig_one_angle(Averaged_sample)...
            + (asin(angle(uv(kk(1),kk(1)))/(2*pi*baseline)))/samples; %Estimated ground angle
    est_sig_one_mag(Averaged_sample) = est_sig_one_mag(Averaged_sample) + abs(uv(kk(1),kk(1)))/samples;     
 
    end
end
 
%% Ploting Results
 
figure(1);
plot(SNR,est_sig_one_angle*180/pi,'bo'); %Estimated ground angle
hold on;
plot(SNR,ones(1,Averaged_samples)*signal_one_angle*180/pi,'b+'); %actual Ground Phase
% plot(SNR,est_sig_two_angle*180/pi,'ro'); %Estimated Vegitation angle
% plot(SNR,ones(1,Noise_samples)*signal_two_angle*180/pi,'r.'); %actual Vegitation phase
title('Direction of Arrival Resolution 1 Delx');
xlabel('SNR dB');ylabel('Angle from referance plane (Degrees)');
legend('Signal One Estimated','Signal One Actual','Location','northeast')
 
axis([-10 20 0 20]);
hold off;
 
figure(2)
title('DOA Coherance 4 Delx');
xlabel('SNR (dB)');ylabel('Maginitude')
hold on;
plot(SNR,est_sig_one_mag,'ro');
legend('Signal Coherance','Location','east');
axis([-10 20 0 1]);
hold off;
