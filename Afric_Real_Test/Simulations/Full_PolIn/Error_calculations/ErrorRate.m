%% Afric_Sim_test_espirit
%% Afric_Sim_test_espirit
clc;clear;
%% Initializations

alpha = 1; %ground weighting factor
pol_signal_one = [1;1;1]/sqrt(3); 
signal_one_offset = 30*pi/180; % ground interferomitry offset
samples = 10;
Window = 100;    %size of window
Averaged_samples = 10000;

signal_one_phase_est = zeros(1,Averaged_samples);
signal_one_mag_est = zeros(1,Averaged_samples);
SNR = zeros(1,Averaged_samples);
%% Matrix Calculations
% Implimenting a window. Esprit and SR techniques
Noise = linspace(.1,10,Averaged_samples);
for Averaged_sample = 1:Averaged_samples;
    
    SNR(Averaged_sample) = 10*log10(1/(Noise(Averaged_sample))^2);
    for n = 1:samples
       
        signal_one =  pol_signal_one*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
        
        s1 = alpha*signal_one;
        s2 = alpha*exp(1i*signal_one_offset)*signal_one ;
        
        s1_Noise = s1 + Noise(Averaged_sample)*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
        s2_Noise = s2 + Noise(Averaged_sample)*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
                                    
        p1 = [s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        p2 = [s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        R1 = p1*p1'/Window;
        R2 = p1*p2'/Window;
        
        A = pinv(R1)*R2;
        
        [u,c] = eig(A);
       
        sort_pol = abs(pol_signal_one'*u);
        [~,kk] = sort(sort_pol,'descend');
        signal_one_phase_est(Averaged_sample) = signal_one_phase_est(Averaged_sample) + angle(c(kk(1),kk(1)))/samples;
        signal_one_mag_est(Averaged_sample) = signal_one_mag_est(Averaged_sample) + abs(c(kk(1),kk(1)))/samples;

    end     
end
%% Plotting Results
signal_one_error = signal_one_phase_est*180/pi - -1*ones(1,Averaged_samples)*signal_one_offset*180/pi;

figure(1);
plot(SNR,signal_one_error,'b')
axis([-20 20 -70 120]);

figure(2)
plot(SNR,-1*signal_one_phase_est*180/pi,'b');
axis([-20 20 -100 100]);

figure(3)
plot(SNR,signal_one_mag_est)