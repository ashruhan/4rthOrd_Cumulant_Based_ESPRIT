%% Afric_Sim_test_ESPRIT
clc;clear;
%% Initializations
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor
eye_optimal = 0.3737;

Pol_ground = [1;-1;0]./sqrt(2);
Pol_vegitation = [1;1;1]./sqrt(3);
G_O = 30;
ground_offset = G_O*pi/180; % ground interferomitry offset
V_O = 50;
vegitation_offset = V_O*pi/180;    % veg interferomitry offset

Averaged_samples = 10;
Window_optimal = 81;    %size of window
SNR_samples = 30;
second_order = 3;

ground_phase_4=zeros(SNR_samples,1);
ground_mag_4=zeros(SNR_samples,1);
ground_phase_2=zeros(SNR_samples,1);
ground_mag_2=zeros(SNR_samples,1);

vegitation_phase_4=zeros(SNR_samples,1);
vegitation_mag_4=zeros(SNR_samples,1);
vegitation_phase_2=zeros(SNR_samples,1);
vegitation_mag_2=zeros(SNR_samples,1);

forth_order = 6;
sort_ground_4 = zeros(2,1);
sort_vegetation_4 = zeros(2,1);

SNR = zeros(1,SNR_samples);
%% Matrix Construction
for SNR_sample = (1:SNR_samples);
    
    SNR(SNR_sample)=SNR_sample-10;
    Noise = (10^(-SNR(SNR_sample)/20))/sqrt(3);
    
    for unusedvariable = 1:Averaged_samples
        
        g =  Pol_ground*(sqrt(-2*log(1-rand(1,Window_optimal))).*exp(1i*2*pi*rand(1,Window_optimal)));
        v =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Window_optimal))).*exp(1i*2*pi*rand(1,Window_optimal)));
        
        s1 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(second_order,Window_optimal))).*exp(1i*2*pi*rand(second_order,Window_optimal));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(second_order,Window_optimal))).*exp(1i*2*pi*rand(second_order,Window_optimal));
        
        
        %% Forming the Six Arrays
        S1_2 = [s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        S2_2 = [s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        R11_norm_2 = S1_2*S1_2.'/Window_optimal;
        E_h1_h1 = R11_norm_2(1,1); E_h1_v1 = R11_norm_2(1,2); E_h1_x1 = R11_norm_2(1,3);
        E_v1_v1 = R11_norm_2(2,2); E_v1_x1 = R11_norm_2(2,3);
        E_x1_x1 = R11_norm_2(3,3);
        
        R11_conj_2 = S1_2*S1_2'/Window_optimal;
        E_h1_h1c = R11_conj_2(1,1); E_h1_v1c = R11_conj_2(1,2); E_h1_x1c = R11_conj_2(1,3);
        E_v1_h1c = R11_conj_2(2,1); E_v1_v1c = R11_conj_2(2,2); E_v1_x1c = R11_conj_2(2,3);
        E_x1_h1c = R11_conj_2(3,1); E_x1_v1c = R11_conj_2(3,2); E_x1_x1c = R11_conj_2(3,3);
        
        R22_norm_2 = S2_2*S2_2.'/Window_optimal;
        E_h2_h2 = R22_norm_2(1,1); E_h2_v2 = R22_norm_2(1,2); E_h2_x2 = R22_norm_2(1,3);
        E_v2_v2 = R22_norm_2(2,2); E_v2_x2 = R22_norm_2(2,3);
        E_x2_x2 = R22_norm_2(3,3);
        
        R12_conj_2 = S1_2*S2_2'/Window_optimal;
        E_h1_h2c = R12_conj_2(1,1); E_h1_v2c = R12_conj_2(1,2); E_h1_x2c = R12_conj_2(1,3);
        E_v1_h2c = R12_conj_2(2,1); E_v1_v2c = R12_conj_2(2,2); E_v1_x2c = R12_conj_2(2,3);
        E_x1_h2c = R12_conj_2(3,1); E_x1_v2c = R12_conj_2(3,2); E_x1_x2c = R12_conj_2(3,3);
        
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
        
        R1_4 = S1_4*S1_4'/Window_optimal;
        R2_4 = S1_4*S2_4'/Window_optimal;
        
        %% Forming Cumulant One Matrix
        Cumulant_1(1,1) = R1_4(1,1) - E_h1_h1*E_h1_h1 - E_h1_h1c*E_h1_h1c -  E_h1_h1c*E_h1_h1c;
        Cumulant_1(1,2) = R1_4(1,2) - E_h1_h1*E_v1_v1 - E_h1_v1c*E_h1_v1c -  E_h1_v1c*E_h1_v1c;
        Cumulant_1(1,3) = R1_4(1,3) - E_h1_h1*E_x1_x1 - E_h1_x1c*E_h1_x1c -  E_h1_x1c*E_h1_x1c;
        Cumulant_1(1,4) = R1_4(1,4) - E_h1_h1*E_h1_v1 - E_h1_h1c*E_h1_v1c -  E_h1_v1c*E_h1_h1c;
        Cumulant_1(1,5) = R1_4(1,5) - E_h1_h1*E_h1_x1 - E_h1_h1c*E_h1_x1c -  E_h1_x1c*E_h1_h1c;
        Cumulant_1(1,6) = R1_4(1,6) - E_h1_h1*E_v1_x1 - E_h1_v1c*E_h1_x1c -  E_h1_x1c*E_h1_v1c;
        
        Cumulant_1(2,1) = R1_4(2,1) - E_v1_v1*E_h1_h1 - E_v1_h1c*E_v1_h1c -  E_v1_h1c*E_v1_h1c;
        Cumulant_1(2,2) = R1_4(2,2) - E_v1_v1*E_v1_v1 - E_v1_v1c*E_v1_v1c -  E_v1_v1c*E_v1_v1c;
        Cumulant_1(2,3) = R1_4(2,3) - E_v1_v1*E_x1_x1 - E_v1_x1c*E_v1_x1c -  E_v1_x1c*E_v1_x1c;
        Cumulant_1(2,4) = R1_4(2,4) - E_v1_v1*E_h1_v1 - E_v1_h1c*E_v1_v1c -  E_v1_v1c*E_v1_h1c;
        Cumulant_1(2,5) = R1_4(2,5) - E_v1_v1*E_h1_x1 - E_v1_h1c*E_v1_x1c -  E_v1_x1c*E_v1_h1c;
        Cumulant_1(2,6) = R1_4(2,6) - E_v1_v1*E_v1_x1 - E_v1_v1c*E_v1_x1c -  E_v1_x1c*E_v1_v1c;
        
        Cumulant_1(3,1) = R1_4(3,1) - E_x1_x1*E_h1_h1 - E_x1_h1c*E_x1_h1c -  E_x1_h1c*E_x1_h1c;
        Cumulant_1(3,2) = R1_4(3,2) - E_x1_x1*E_v1_v1 - E_x1_v1c*E_x1_v1c -  E_x1_v1c*E_x1_v1c;
        Cumulant_1(3,3) = R1_4(3,3) - E_x1_x1*E_x1_x1 - E_x1_x1c*E_x1_x1c -  E_x1_x1c*E_x1_x1c;
        Cumulant_1(3,4) = R1_4(3,4) - E_x1_x1*E_h1_v1 - E_x1_h1c*E_x1_v1c -  E_x1_v1c*E_x1_h1c;
        Cumulant_1(3,5) = R1_4(3,5) - E_x1_x1*E_h1_x1 - E_x1_h1c*E_x1_x1c -  E_x1_x1c*E_x1_h1c;
        Cumulant_1(3,6) = R1_4(3,6) - E_x1_x1*E_v1_x1 - E_x1_v1c*E_x1_x1c -  E_x1_x1c*E_x1_v1c;
        
        Cumulant_1(4,1) = R1_4(4,1) - E_h1_v1*E_h1_h1 - E_h1_h1c*E_v1_h1c -  E_h1_h1c*E_v1_h1c;
        Cumulant_1(4,2) = R1_4(4,2) - E_h1_v1*E_v1_v1 - E_h1_v1c*E_v1_v1c -  E_h1_v1c*E_v1_v1c;
        Cumulant_1(4,3) = R1_4(4,3) - E_h1_v1*E_x1_x1 - E_h1_x1c*E_v1_x1c -  E_h1_x1c*E_v1_x1c;
        Cumulant_1(4,4) = R1_4(4,4) - E_h1_v1*E_h1_v1 - E_h1_h1*E_v1_v1c  -  E_h1_v1c*E_v1_h1c;
        Cumulant_1(4,5) = R1_4(4,5) - E_h1_v1*E_h1_x1 - E_h1_h1c*E_v1_x1c -  E_h1_x1c*E_v1_h1c;
        Cumulant_1(4,6) = R1_4(4,6) - E_h1_v1*E_v1_x1 - E_h1_v1c*E_v1_x1c -  E_h1_x1c*E_v1_v1c;
        
        Cumulant_1(5,1) = R1_4(5,1) - E_h1_x1*E_h1_h1 - E_h1_h1c*E_x1_h1c -  E_h1_h1c*E_x1_h1c;
        Cumulant_1(5,2) = R1_4(5,2) - E_h1_x1*E_v1_v1 - E_h1_v1c*E_x1_v1c -  E_h1_v1c*E_x1_v1c;
        Cumulant_1(5,3) = R1_4(5,3) - E_h1_x1*E_x1_x1 - E_h1_x1c*E_x1_x1c -  E_h1_x1c*E_x1_x1c;
        Cumulant_1(5,4) = R1_4(5,4) - E_h1_x1*E_h1_v1 - E_h1_h1c*E_x1_v1c -  E_h1_v1c*E_x1_h1c;
        Cumulant_1(5,5) = R1_4(5,5) - E_h1_x1*E_h1_x1 - E_h1_h1c*E_x1_x1c -  E_h1_x1c*E_x1_h1c;
        Cumulant_1(5,6) = R1_4(5,6) - E_h1_x1*E_v1_x1 - E_h1_v1c*E_x1_x1c -  E_h1_x1c*E_x1_v1c;
        
        Cumulant_1(6,1) = R1_4(6,1) - E_v1_x1*E_h1_h1 - E_v1_h1c*E_x1_h1c -  E_v1_h1c*E_x1_h1c;
        Cumulant_1(6,2) = R1_4(6,2) - E_v1_x1*E_v1_v1 - E_v1_v1c*E_x1_v1c -  E_v1_v1c*E_x1_v1c;
        Cumulant_1(6,3) = R1_4(6,3) - E_v1_x1*E_x1_x1 - E_v1_x1c*E_x1_x1c -  E_v1_x1c*E_x1_x1c;
        Cumulant_1(6,4) = R1_4(6,4) - E_v1_x1*E_h1_v1 - E_v1_h1c*E_x1_v1c -  E_v1_v1c*E_x1_h1c;
        Cumulant_1(6,5) = R1_4(6,5) - E_v1_x1*E_h1_x1 - E_v1_h1c*E_x1_x1c -  E_v1_x1c*E_x1_h1c;
        Cumulant_1(6,6) = R1_4(6,6) - E_v1_x1*E_v1_x1 - E_v1_v1c*E_x1_x1c -  E_v1_x1c*E_x1_v1c;
        
        %% Forming Cumulant Two Matrix
        Cumulant_2(1,1) = R2_4(1,1) - E_h1_h1*E_h2_h2 - E_h1_h2c*E_h1_h2c -  E_h1_h2c*E_h1_h2c;
        Cumulant_2(1,2) = R2_4(1,2) - E_h1_h1*E_v2_v2 - E_h1_v2c*E_h1_v2c -  E_h1_v2c*E_h1_v2c;
        Cumulant_2(1,3) = R2_4(1,3) - E_h1_h1*E_x2_x2 - E_h1_x2c*E_h1_x2c -  E_h1_x2c*E_h1_x2c;
        Cumulant_2(1,4) = R2_4(1,4) - E_h1_h1*E_h2_v2 - E_h1_h2c*E_h1_v2c -  E_h1_v2c*E_h1_h2c;
        Cumulant_2(1,5) = R2_4(1,5) - E_h1_h1*E_h2_x2 - E_h1_h2c*E_h1_x2c -  E_h1_x2c*E_h1_h2c;
        Cumulant_2(1,6) = R2_4(1,6) - E_h1_h1*E_v2_x2 - E_h1_v2c*E_h1_x2c -  E_h1_x2c*E_h1_v2c;
        
        Cumulant_2(2,1) = R2_4(2,1) - E_v1_v1*E_h2_h2 - E_v1_h2c*E_v1_h2c -  E_v1_h2c*E_v1_h2c;
        Cumulant_2(2,2) = R2_4(2,2) - E_v1_v1*E_v2_v2 - E_v1_v2c*E_v1_v2c -  E_v1_v2c*E_v1_v2c;
        Cumulant_2(2,3) = R2_4(2,3) - E_v1_v1*E_x2_x2 - E_v1_x2c*E_v1_x2c -  E_v1_x2c*E_v1_x2c;
        Cumulant_2(2,4) = R2_4(2,4) - E_v1_v1*E_h2_v2 - E_v1_h2c*E_v1_v2c -  E_v1_v2c*E_v1_h2c;
        Cumulant_2(2,5) = R2_4(2,5) - E_v1_v1*E_h2_x2 - E_v1_h2c*E_v1_x2c -  E_v1_x2c*E_v1_h2c;
        Cumulant_2(2,6) = R2_4(2,6) - E_v1_v1*E_v2_x2 - E_v1_v2c*E_v1_x2c -  E_v1_x2c*E_v1_v2c;
        
        Cumulant_2(3,1) = R2_4(3,1) - E_x1_x1*E_h2_h2 - E_x1_h2c*E_x1_h2c -  E_x1_h2c*E_x1_h2c;
        Cumulant_2(3,2) = R2_4(3,2) - E_x1_x1*E_v2_v2 - E_x1_v2c*E_x1_v2c -  E_x1_v2c*E_x1_v2c;
        Cumulant_2(3,3) = R2_4(3,3) - E_x1_x1*E_x2_x2 - E_x1_x2c*E_x1_x2c -  E_x1_x2c*E_x1_x2c;
        Cumulant_2(3,4) = R2_4(3,4) - E_x1_x1*E_h2_v2 - E_x1_h2c*E_x1_v2c -  E_x1_v2c*E_x1_h2c;
        Cumulant_2(3,5) = R2_4(3,5) - E_x1_x1*E_h2_x2 - E_x1_h2c*E_x1_x2c -  E_x1_x2c*E_x1_h2c;
        Cumulant_2(3,6) = R2_4(3,6) - E_x1_x1*E_v2_x2 - E_x1_v2c*E_x1_x2c -  E_x1_x2c*E_x1_v2c;
        
        Cumulant_2(4,1) = R2_4(4,1) - E_h1_v1*E_h2_h2 - E_h1_h2c*E_v1_h2c -  E_h1_h2c*E_v1_h2c;
        Cumulant_2(4,2) = R2_4(4,2) - E_h1_v1*E_v2_v2 - E_h1_v2c*E_v1_v2c -  E_h1_v2c*E_v1_v2c;
        Cumulant_2(4,3) = R2_4(4,3) - E_h1_v1*E_x2_x2 - E_h1_x2c*E_v1_x2c -  E_h1_x2c*E_v1_x2c;
        Cumulant_2(4,4) = R2_4(4,4) - E_h1_v1*E_h2_v2 - E_h1_h2c*E_v1_v2c -  E_h1_v2c*E_v1_h2c;
        Cumulant_2(4,5) = R2_4(4,5) - E_h1_v1*E_h2_x2 - E_h1_h2c*E_v1_x2c -  E_h1_x2c*E_v1_h2c;
        Cumulant_2(4,6) = R2_4(4,6) - E_h1_v1*E_v2_x2 - E_h1_v2c*E_v1_x2c -  E_h1_x2c*E_v1_v2c;
        
        Cumulant_2(5,1) = R2_4(5,1) - E_h1_x1*E_h2_h2 - E_h1_h2c*E_x1_h2c -  E_h1_h2c*E_x1_h2c;
        Cumulant_2(5,2) = R2_4(5,2) - E_h1_x1*E_v2_v2 - E_h1_v2c*E_x1_v2c -  E_h1_v2c*E_x1_v2c;
        Cumulant_2(5,3) = R2_4(5,3) - E_h1_x1*E_x2_x2 - E_h1_x2c*E_x1_x2c -  E_h1_x2c*E_x1_x2c;
        Cumulant_2(5,4) = R2_4(5,4) - E_h1_x1*E_h2_v2 - E_h1_h2c*E_x1_v2c -  E_h1_v2c*E_x1_h2c;
        Cumulant_2(5,5) = R2_4(5,5) - E_h1_x1*E_h2_x2 - E_h1_h2c*E_x1_x2c -  E_h1_x2c*E_x1_h2c;
        Cumulant_2(5,6) = R2_4(5,6) - E_h1_x1*E_v2_x2 - E_h1_v2c*E_x1_x2c -  E_h1_x2c*E_x1_v2c;
        
        Cumulant_2(6,1) = R2_4(6,1) - E_v1_x1*E_h2_h2 - E_v1_h2c*E_x1_h2c -  E_v1_h2c*E_x1_h2c;
        Cumulant_2(6,2) = R2_4(6,2) - E_v1_x1*E_v2_v2 - E_v1_v2c*E_x1_v2c -  E_v1_v2c*E_x1_v2c;
        Cumulant_2(6,3) = R2_4(6,3) - E_v1_x1*E_x2_x2 - E_v1_x2c*E_x1_x2c -  E_v1_x2c*E_x1_x2c;
        Cumulant_2(6,4) = R2_4(6,4) - E_v1_x1*E_h2_v2 - E_v1_h2c*E_x1_v2c -  E_v1_v2c*E_x1_h2c;
        Cumulant_2(6,5) = R2_4(6,5) - E_v1_x1*E_h2_x2 - E_v1_h2c*E_x1_x2c -  E_v1_x2c*E_x1_h2c;
        Cumulant_2(6,6) = R2_4(6,6) - E_v1_x1*E_v2_x2 - E_v1_v2c*E_x1_x2c -  E_v1_x2c*E_x1_v2c;
        
        %% ESPRIT Algorithm
        [eigenvec_4,eigenval_4] = eig(pinv(Cumulant_1)*Cumulant_2);
        
        for i = 1:forth_order
            sort_ground_4(i) = (abs(eigenval_4(i,i)))^2....
                *(abs(eigenvec_4(1,i))^2....
                + abs(eigenvec_4(2,i))^2....
                + abs(eigenvec_4(4,i))^2);
            
            sort_vegetation_4(i) = (abs(eigenval_4(i,i)))^2....
                *(abs(eigenvec_4(3,i))^2....
                + abs(eigenvec_4(5,i))^2....
                + abs(eigenvec_4(6,i))^2);
        end
        
        [~,srt_g_4] = sort(sort_ground_4,'descend');
        [~,srt_v_4] = sort(sort_vegetation_4,'descend');
        
        ground_phase_4(SNR_sample) = ground_phase_4(SNR_sample)....
            + 0.5*angle(eigenval_4(srt_g_4(1),srt_g_4(1)))/Averaged_samples;
        
        vegitation_phase_4(SNR_sample) = vegitation_phase_4(SNR_sample)....
            + 0.5*angle(eigenval_4(srt_v_4(1),srt_v_4(1)))/Averaged_samples;
    end
end
%% Plotting Results
figure(1);
title('2nd and 4rth Order Modified ESPRIT Interferometric Phases');
xlabel('SNR dB');ylabel('Int Phase (Degrees)');
hold on;
plot(SNR,ground_phase_4*180/pi,'bx');
plot(SNR,vegitation_phase_4*180/pi,'gx');
plot(SNR,-V_O*ones(1,SNR_samples),'g');
plot(SNR,-G_O*ones(1,SNR_samples),'b');
legend('4rth Order Ground','2nd Order Ground','4rth Order Vegitation','2nd Order Vegitaion','Location','northeast')
hold off