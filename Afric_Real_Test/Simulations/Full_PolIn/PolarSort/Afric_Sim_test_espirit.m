%% Afric_Sim_test_espirit
clc;clear;
%% Initializations
% Setting up the enviornment
g_weight = 1; %ground weighting factor
v_weight = 1;   %veg weighting factor
n_weight = 1*10^(-3/2);% Added Noise to the System

Pol_ground = [1;1;0]; %Multivatiant Ground
Pol_vegitation = [1;0;1]; %Multivariant Vegitation

ground_offset = pi/4; % Ground Interferomitry offset
vegitation_offset = pi/3;    % Vegitation Interferomitry offset

Averaging_loop_size = 50;
Signal_samples = 100;    %size of Ensamble Average Window

Noise_samples = 50;
est_ground_angle = zeros(Noise_samples,1); est_vegitation_angle = zeros(Noise_samples,1);
%% Matrix Calculations
% Implimenting a window. Esprit and SR techniques
for Averaged_sample = 1:Noise_samples;
    Noise = Averaged_sample*n_weight;
    for unusedVariable = 1:Averaging_loop_size;
        
        g =  Pol_ground*(sqrt(-2*log(1-rand(1,Signal_samples))).*exp(1i*2*pi*rand(1,Signal_samples)));
        v =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Signal_samples))).*exp(1i*2*pi*rand(1,Signal_samples)));
        
        signal_1 = g_weight*g + v_weight*v;
        signal_2 = g_weight*exp(1i*ground_offset)*g + g_weight*exp(1i*vegitation_offset)*v;
        
        AddedNoise_1 = Noise*sqrt(-2*log(1-rand(3,Signal_samples))).*exp(1i*2*pi*rand(3,Signal_samples));
        AddedNoise_2 = Noise*sqrt(-2*log(1-rand(3,Signal_samples))).*exp(1i*2*pi*rand(3,Signal_samples));

        S1 = signal_1 + AddedNoise_1;
        S2 = signal_2 + AddedNoise_2;
        
        R1 = S1*S1';
        R2 = S1*S2';
        
        A = pinv(R1)*R2;
        [u,c] = eig(A);
        
        sg = abs(Pol_ground'*u);
        [~,kk] = sort(sg,'descend');
        est_ground_angle(Averaged_sample) = est_ground_angle(Averaged_sample) + angle(c(kk(1),kk(1)))/Averaging_loop_size;
        
        sv = abs(Pol_vegitation'*u);
        [~,kk] = sort(sv,'descend');
        est_vegitation_angle(Averaged_sample) = est_vegitation_angle(Averaged_sample) + angle(c(kk(1),kk(1)))/Averaging_loop_size;
        
    end
end
temp = 1:Noise_samples;
n_snr = 10*log10(1./(temp*n_weight).^2)';
%% Ploting Results
g_error = est_ground_angle*180/pi - -1*ones(Noise_samples,1)*ground_offset*180/pi;
v_error = est_vegitation_angle*180/pi - -1*ones(Noise_samples,1)*vegitation_offset*180/pi;
figure(1);title('Ground in blue Veg in red');
hold on;
plot(n_snr,g_error,'b')
plot(n_snr,v_error,'r')
xlabel('SNR dB');ylabel('error in degrees')
hold off;

figure(2);title('Ground and Vegitation angle in degrees');
plot(n_snr,est_ground_angle*180/pi,'bo'); %Estimated ground angle
hold on;
plot(n_snr,-1*ones(Noise_samples,1)*ground_offset*180/pi,'b.'); %actual Ground Phase
plot(n_snr,est_vegitation_angle*180/pi,'ro'); %Estimated Vegitation angle
plot(n_snr,-1*ones(Noise_samples,1)*vegitation_offset*180/pi,'r.'); %actual Vegitation phase
hold off;