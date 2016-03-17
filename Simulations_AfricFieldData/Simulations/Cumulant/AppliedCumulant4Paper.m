%% Afric_Sim_test_ESPRIT
clc;clear;
%% Initializations
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor

Pol_ground = [1;1;0]./sqrt(2);
Pol_vegitation = [1;1;1]./sqrt(3);

eye_4 = -0.099;
eye_2 = 0;
G_O = 20;
ground_offset = G_O*pi/180; % ground interferomitry offset
V_O = 60;
vegitation_offset = V_O*pi/180;    % veg interferomitry offset

Averaged_samples = 20;
Window_optimal = 81;    %size of window
SNR_samples = 30;

ground_angle_4 = zeros(SNR_samples,1);
ground_angle_2 = zeros(SNR_samples,1);

vegitation_angle_4 = zeros(SNR_samples,1);
vegitation_angle_2 = zeros(SNR_samples,1);

ground_abs_4 = zeros(SNR_samples,1);
ground_abs_2 = zeros(SNR_samples,1);

vegitation_abs_4 = zeros(SNR_samples,1);
vegitation_abs_2 = zeros(SNR_samples,1);

ground_angle_rmse_4 = zeros(SNR_samples,1);
vegitation_angle_rmse_4 = zeros(SNR_samples,1);

ground_angle_rmse_2 = zeros(SNR_samples,1);
vegitation_angle_rmse_2 = zeros(SNR_samples,1);

SNR = zeros(1,SNR_samples);
%% Matrix Construction
for SNR_sample = fliplr(1:SNR_samples);
    
    SNR(SNR_sample)=SNR_sample-15;
    Noise = (10^(-SNR(SNR_sample)/20))/sqrt(3);
    
    for unusedvariable = 1:Averaged_samples
        
        %% Matrix Construction
        
        g =  Pol_ground*exp(1i*2*pi*rand(1,Window_optimal));
        v =  Pol_vegitation*exp(1i*2*pi*rand(1,Window_optimal));
        
        s1 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(3,Window_optimal))).*exp(1i*2*pi*rand(3,Window_optimal));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,Window_optimal))).*exp(1i*2*pi*rand(3,Window_optimal));
        
        %% Second Order Stats
        
        S1_2 = [s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        S2_2 = [s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        R1_2 = S1_2*S1_2'/Window_optimal;
        R2_2 = S1_2*S2_2'/Window_optimal;
        
        [eigenvec_2,eigenval_2] =eig((pinv(R1_2 + eye_2*eye(3)))*R2_2,'nobalance');
        
        [~,srt_2]=sort(abs(diag(eigenval_2)),'descend');
        
        Leig_copol = (abs(eigenval_2(srt_2(1),srt_2(1)))^2)....
            *abs(eigenvec_2(1,srt_2(1)))^2....
            + abs(eigenvec_2(2,srt_2(1)))^2;
        
        SLeig_copol = (abs(eigenval_2(srt_2(1),srt_2(1)))^2)....
            *abs(eigenvec_2(1,srt_2(2)))^2....
            + abs(eigenvec_2(2,srt_2(2)))^2;
        
        if (Leig_copol >= SLeig_copol)
            
            ground_angle_2(SNR_sample,1) = ground_angle_2(SNR_sample,1)....
                + angle(eigenval_2(srt_2(1),srt_2(1)))/Averaged_samples;
            
            ground_angle_rmse_2(SNR_sample) =  ground_angle_rmse_2(SNR_sample)...
                + ((ground_offset + angle(eigenval_2(srt_2(1),srt_2(1))))^2)/Averaged_samples;
            
            
            vegitation_angle_2(SNR_sample,1) = vegitation_angle_2(SNR_sample,1)....
                + angle(eigenval_2(srt_2(2),srt_2(2)))/Averaged_samples;
            
            vegitation_angle_rmse_2(SNR_sample) =  vegitation_angle_rmse_2(SNR_sample)...
                + ((vegitation_offset + angle(eigenval_2(srt_2(2),srt_2(2))))^2)/Averaged_samples;
            
            
            
            ground_abs_2(SNR_sample,1) = ground_abs_2(SNR_sample,1)....
                + (abs(eigenval_2(srt_2(1),srt_2(1))))/Averaged_samples;
            
            vegitation_abs_2(SNR_sample,1) = vegitation_abs_2(SNR_sample,1)....
                + (abs(eigenval_2(srt_2(2),srt_2(2))))/Averaged_samples;
        else
            
            ground_angle_2(SNR_sample,1) = ground_angle_2(SNR_sample,1)....
                + angle(eigenval_2(srt_2(2),srt_2(2)))/Averaged_samples;
            
            ground_angle_rmse_2(SNR_sample) =  ground_angle_rmse_2(SNR_sample)...
                + ((ground_offset + angle(eigenval_2(srt_2(2),srt_2(2))))^2)/Averaged_samples;
            
            
            vegitation_angle_2(SNR_sample,1) = vegitation_angle_2(SNR_sample,1)....
                + angle(eigenval_2(srt_2(1),srt_2(1)))/Averaged_samples;
            
            vegitation_angle_rmse_2(SNR_sample) =  vegitation_angle_rmse_2(SNR_sample)...
                + ((vegitation_offset + angle(eigenval_2(srt_2(1),srt_2(1))))^2)/Averaged_samples;
            
            
            ground_abs_2(SNR_sample,1) = ground_abs_2(SNR_sample,1)....
                +(abs(eigenval_2(srt_2(2),srt_2(2))))/Averaged_samples;
            
            vegitation_abs_2(SNR_sample,1) = vegitation_abs_2(SNR_sample,1)....
                + (abs(eigenval_2(srt_2(1),srt_2(1))))/Averaged_samples;
        end
        
        %% Fourth Order Statistics
        [ Cumulant_11, Cumulant_12] = Cumulant( s1_Noise ,s2_Noise,Window_optimal );
        
        [eigenvec_4,eigenval_4] = eig((pinv(Cumulant_11+eye_4*eye(6)))...
            *Cumulant_12,'nobalance');
        [~,srt_4] = sort(abs(diag(eigenval_4)),'descend');
        
        LeigTemp  = (abs(eigenval_4(srt_4(1),srt_4(1))))^2....
            *(abs(eigenvec_4(3,srt_4(1)))^2....
            + abs(eigenvec_4(5,srt_4(1)))^2....
            + abs(eigenvec_4(6,srt_4(1)))^2);
        
        SLeigTemp = (abs(eigenval_4(srt_4(2),srt_4(2))))^2....
            *(abs(eigenvec_4(3,srt_4(2)))^2....
            + abs(eigenvec_4(5,srt_4(2)))^2....
            + abs(eigenvec_4(6,srt_4(2)))^2);
        
        if LeigTemp >= SLeigTemp
            vegitation_angle_4(SNR_sample,1) = vegitation_angle_4(SNR_sample,1)....
                + 0.5*angle(eigenval_4(srt_4(1),srt_4(1)))/Averaged_samples;
            
            vegitation_angle_rmse_4(SNR_sample) =  vegitation_angle_rmse_4(SNR_sample)...
                + ((vegitation_offset - abs(0.5*angle(eigenval_4(srt_4(1),srt_4(1)))))^2)/Averaged_samples;
            
            ground_angle_4(SNR_sample,1) = ground_angle_4(SNR_sample,1)....
                + 0.5*angle(eigenval_4(srt_4(2),srt_4(2)))/Averaged_samples;
            
            ground_angle_rmse_4(SNR_sample) =  ground_angle_rmse_4(SNR_sample)...
                + ((ground_offset - abs(0.5*angle(eigenval_4(srt_4(2),srt_4(2)))))^2)/Averaged_samples;
            
            
            vegitation_abs_4(SNR_sample,1) = vegitation_abs_4(SNR_sample,1)....
                + sqrt(abs(eigenval_4(srt_4(1),srt_4(1))))/Averaged_samples;
            
            ground_abs_4(SNR_sample,1) = ground_abs_4(SNR_sample,1)....
                + sqrt(abs(eigenval_4(srt_4(2),srt_4(2))))/Averaged_samples;
            
        else
            
            vegitation_angle_4(SNR_sample,1) = vegitation_angle_4(SNR_sample,1)....
                + 0.5*angle(eigenval_4(srt_4(2),srt_4(2)))/Averaged_samples;
            
            vegitation_angle_rmse_4(SNR_sample) =  vegitation_angle_rmse_4(SNR_sample)...
                + ((vegitation_offset - abs(0.5*angle(eigenval_4(srt_4(2),srt_4(2)))))^2)/Averaged_samples;
            
            
            ground_angle_4(SNR_sample,1) = ground_angle_4(SNR_sample,1)....
                + 0.5*angle(eigenval_4(srt_4(1),srt_4(1)))/Averaged_samples;
            
            ground_angle_rmse_4(SNR_sample) =  ground_angle_rmse_4(SNR_sample)...
                + ((ground_offset - abs(0.5*angle(eigenval_4(srt_4(1),srt_4(1)))))^2)/Averaged_samples;
            
            
            vegitation_abs_4(SNR_sample,1) = vegitation_abs_4(SNR_sample,1)....
                + sqrt(abs(eigenval_4(srt_4(2),srt_4(2))))/Averaged_samples;
            
            ground_abs_4(SNR_sample,1) = ground_abs_4(SNR_sample,1)....
                + sqrt(abs(eigenval_4(srt_4(1),srt_4(1))))/Averaged_samples;
        end
    end
end
ground_angle_rmse_2 = sqrt(ground_angle_rmse_2)*180/pi;
vegitation_angle_rmse_2 = sqrt(vegitation_angle_rmse_2)*180/pi;
ground_angle_rmse_4 = sqrt(ground_angle_rmse_4)*180/pi;
vegitation_angle_rmse_4 = sqrt(vegitation_angle_rmse_4)*180/pi;
%% Plotting Results

figure(1);
title('4rth Order ESPRITInterferometric Phases');
xlabel('SNR (dB)');ylabel('Interferometric Phase (Degrees)');
hold on;
plot(SNR,(-ground_angle_4)*180/pi,'bx');
plot(SNR,(-vegitation_angle_4)*180/pi,'gx');
plot(SNR,G_O*ones(1,SNR_samples),'b');
plot(SNR,V_O*ones(1,SNR_samples),'g');
axis([-15,15,G_O-5,V_O+5])
legend('4rth Order Ground Estimate','4rth Order Vegetation Estimate',...
    '4rth Order Ground Actual','4rth Order Vegetation Actual','Location','east')
hold off

figure(2);
title('2nd Order ESPRIT Interferometric Phases');
xlabel('SNR (dB)');ylabel('Interferometric Phase (Degrees)');
hold on;
plot(SNR,(-ground_angle_2)*180/pi,'bo');
plot(SNR,(-vegitation_angle_2)*180/pi,'go');
plot(SNR,G_O*ones(1,SNR_samples),'b');
plot(SNR,V_O*ones(1,SNR_samples),'g');
axis([-15,15,G_O-5,V_O+5])
legend('2nd Order Ground Estimate','2nd Order Vegetation Estimate',...
    '2nd Order Ground Actual','2nd Order Vegetation Actual','Location','east')
hold off
%%
figure(3);
title('4rth Order ESPRIT Coherance');
xlabel('SNR (dB)');ylabel('Magnitude');
hold on;
plot(SNR,ground_abs_4,'bx');
plot(SNR,vegitation_abs_4,'gx');
axis([-15,15,0,2])
legend('4rth Order Ground','4rth Order Vegetation','Location','northwest')
hold off
%%
figure(4);
title('2nd Order ESPRIT Coherance');
xlabel('SNR (dB)');ylabel('Magnitude');
hold on;
plot(SNR,ground_abs_2,'bo');
plot(SNR,vegitation_abs_2,'go');
legend('2nd Order Ground','2nd Order Vegetaion','Location','northwest')
hold off

figure(5)
title('4rth Order ESPRIT RMS Error vs SNR');
xlabel('SNR (dB)');ylabel('RMS Error (degrees)');
hold on;
plot(SNR,ground_angle_rmse_4,'bx')
plot(SNR,vegitation_angle_rmse_4,'gx')
legend('4rth Order Ground','4rth Order Vegetation','Location','northwest')
hold off

figure(6)
title('2nd Order ESPRIT RMS Error vs SNR');
xlabel('SNR (dB)');ylabel('RMS Error (degrees)');
hold on;
plot(SNR,ground_angle_rmse_2,'bo')
plot(SNR,vegitation_angle_rmse_2,'go')
legend('2nd Order Ground','2nd Order Vegetaion','Location','northwest')
hold off
%%
figure(7)
title('4rth Order ESPRIT RMS Error vs Coherence');
xlabel('RMS Error (degrees)');ylabel('Magnitude');
hold on;
plot(ground_angle_rmse_4,ground_abs_4,'bx')
plot(vegitation_angle_rmse_4,vegitation_abs_4,'gx')
legend('4rth Order Ground','4rth Order Vegetaion','Location','northwest');
legend('4rth Order Ground','4rth Order Vegetation','Location','northwest');
hold off
%%
figure(8)
title('2nd Order ESPRIT RMS Error vs Coherence');
xlabel('RMS Error (degrees)');ylabel('Magnitude');
hold on;
plot(ground_angle_rmse_2,ground_abs_2,'bx')
plot(vegitation_angle_rmse_2,vegitation_abs_2,'gx')
legend('2nd Order Ground','2nd Order Vegetaion','Location','northwest')
hold off

