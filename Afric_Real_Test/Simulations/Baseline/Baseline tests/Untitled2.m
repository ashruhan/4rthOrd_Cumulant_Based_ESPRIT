%% Espirit_baseline
clear;clc;
%% Initializations
g_weight = 1;
v_weight = 1;
n_weight = 1*10^(-3/2);

g_angle = -pi/24;
v_angle = pi/16;

delx = 0.5; % Delta x of baseline
baseline = delx; %Physical baseline of two sample arrays

g_phase = sin(g_angle)*2*pi*baseline; %Used in the sorting routine
v_phase = sin(v_angle)*2*pi*baseline; %Used in the sorting routine

array_size = 6; % the amount of ground samples
Averaging_loop_size = 100;
Signal_samples = 20; %Averaged_Signal_samples in matrix

Noise_samples = 20;
est_ground_angle = zeros(1,Noise_samples);
est_vegitation_angle = zeros(1,Noise_samples);

S1Pos1 = zeros(Signal_samples,Signal_samples);
S1Pos2 = zeros(Signal_samples,Signal_samples);
S2Pos1 = zeros(Signal_samples,Signal_samples);
S2Pos2 = zeros(Signal_samples,Signal_samples);
%% Start of Algorithm
for Noise_Increment = 1:Noise_samples;
    Noise_W = Noise_Increment*n_weight;
    for unusedVariable = 1:Averaging_loop_size;
        %% Setting up Random Phi that differ across samples
        phi_one = 1i*2*pi*rand(1,Signal_samples);
        phi_two = 1i*2*pi*rand(1,Signal_samples);
        pos = [1:array_size]'; %#ok<NBRAK>
        distance = 2.*pi.*pos.*delx;  % element separation of 1/2 wavelength
        
        %% Element Array size with Signal
        
        X = g_weight.*(exp(-1i*distance*sin(g_angle))*exp(phi_one))... %Ground
            + v_weight.*(exp(-1i*distance*sin(v_angle))*exp(phi_two))...  %Vegetaion
            + Noise_W.*sqrt(-2.*log(1-rand(array_size,Signal_samples))).*exp(1i*2*pi*rand(array_size,Signal_samples));  %Noise
        
        %% Setting up the Shifting Array  
         % Setting up Baselines
         
        % Remember that if you change baseline to change the Initial
        % Varialble to reflect changes
         
        for index = 1:Signal_samples
            
            if index == 1
                 S1Pos1(1,:) = X(1,:); S1Pos2(1,:) = X(2,:); S2Pos1(1,:) = X(5,:); S2Pos2(1,:) = X(6,:); %Baseline = 4delx
%                 S1Pos1(1,:) = X(2,:); S1Pos2(1,:) = X(3,:); S2Pos1(1,:) = X(4,:); S2Pos2(1,:) = X(5,:); %Baseline = 2delx
%                 S1Pos1(1,:) = X(3,:); S1Pos2(1,:) = X(4,:); S2Pos1(1,:) = X(4,:); S2Pos2(1,:) = X(5,:); %Baseline = delx
            else
                S1Pos1(index,:) = circshift(S1Pos1(index-1,:),[0,1]);
                S1Pos2(index,:) = circshift(S1Pos2(index-1,:),[0,1]);
                S2Pos1(index,:) = circshift(S2Pos1(index-1,:),[0,1]);
                S2Pos2(index,:) = circshift(S2Pos2(index-1,:),[0,1]);
            end
            
        end        
        %% Espirit Algorithm
                
        S1 = [S1Pos1; S1Pos2];
        S2 = [S2Pos1; S2Pos2];
        
        R1 = S1*S1';
        R2 = S1*S2';
        A =  pinv(R1)*R2;
        
        %% Averaging over a matix then averaging that signal in a for loop proportional to the matix rows
        
        [u,uv] = eig(A);
        
        [~,kk]=sort(abs(angle(diag(uv))-g_phase),'ascend');
        est_ground_angle(Noise_Increment) = est_ground_angle(Noise_Increment)...
            + (asin(angle(uv(kk(1),kk(1)))/(2*pi*baseline)))/Averaging_loop_size; %Estimated ground angle
        
        [~,kk]=sort(abs(angle(diag(uv))-v_phase),'ascend');
        est_vegitation_angle(Noise_Increment) = est_vegitation_angle(Noise_Increment)...
            + (asin(angle(uv(kk(1),kk(1)))/(2*pi*baseline)))/Averaging_loop_size; %Estimated Vegitation angle
    end
end

%% Noise Power
temp = 1:Noise_samples;
n_snr = 10*log10(1./(temp*n_weight).^2);
%% Ploting Results
g_error = est_ground_angle*180/pi - ones(1,Noise_samples)*g_angle*180/pi;
v_error = est_vegitation_angle*180/pi - ones(1,Noise_samples)*v_angle*180/pi;
figure(3);title('Ground in blue Veg in red');
hold on;
plot(n_snr,g_error,'b')
plot(n_snr,v_error,'r')
xlabel('SNR dB');ylabel('error in degrees')
hold off;

figure(4);title('Ground in blue Veg in red');
plot(n_snr,est_ground_angle*180/pi,'bo'); %Estimated ground angle
hold on;
plot(n_snr,ones(1,Noise_samples)*g_angle*180/pi,'b.'); %actual Ground Phase
plot(n_snr,est_vegitation_angle*180/pi,'ro'); %Estimated Vegitation angle
plot(n_snr,ones(1,Noise_samples)*v_angle*180/pi,'r.'); %actual Vegitation phase
xlabel('SNR dB');ylabel('error in degrees')
hold off;