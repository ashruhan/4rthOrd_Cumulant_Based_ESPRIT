%% Afric_Sim_test_ESPRIT
clc;clear;
%% Initializations
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor

Pol_ground = [1;-1;0]./sqrt(2);
Pol_vegitation = [1;1;1]./sqrt(3);

G_O = 20;
ground_offset = G_O*pi/180; % ground interferomitry offset
V_O = 60;
vegitation_offset = V_O*pi/180;    % veg interferomitry offset

Averaged_samples = 100;
Window_optimal = 81;    %size of window
SNR_samples = 30;

ground_angle_4 = zeros(SNR_samples,1);
vegitation_angle_4 = zeros(SNR_samples,1);
ground_abs_4 = zeros(SNR_samples,1);
vegitation_abs_4 = zeros(SNR_samples,1);
ground_angle_rmse_4 = zeros(SNR_samples,1);
vegitation_angle_rmse_4 = zeros(SNR_samples,1);

SNR = zeros(1,SNR_samples);
%% Matrix Construction
for SNR_sample = fliplr(1:SNR_samples);
    
    SNR(SNR_sample)=SNR_sample-15;
    Noise = (10^(-SNR(SNR_sample)/20))/sqrt(3);
    
    for unusedvariable = 1:Averaged_samples
        
        %% Signal Construction

        g =  Pol_ground*(sqrt(-2*log(1-rand(1,Window_optimal))).*exp(1i*2*pi*rand(1,Window_optimal)));
        v =  Pol_vegitation*(sqrt(-2*log(1-rand(1,Window_optimal))).*exp(1i*2*pi*rand(1,Window_optimal)));
        
        s1 = alpha*g + beta*v;
        s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegitation_offset)*v;
        
        s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(3,Window_optimal))).*exp(1i*2*pi*rand(3,Window_optimal));
        s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,Window_optimal))).*exp(1i*2*pi*rand(3,Window_optimal));
        
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
        
        %% Fourth Order Statistics
        clc;
        [~,eigenvalCov11_4] = eig(R1_4,'nobalance');
        
        [~,srtCov_4] = sort(abs(diag(eigenvalCov11_4)),'descend');
        
        eye_4 = mean((diag(eigenvalCov11_4)))/max(abs(diag(eigenvalCov11_4)))^2;
        
        [eigenvec_4,eigenval_4] = eig((pinv(R1_4 + eye_4*eye(6)))*R2_4,'nobalance');
        
        [~,srt_abs_4] = sort(abs(diag(eigenval_4)),'descend');
        
        [~,srt_angle_4] = sort(abs(angle(diag(eigenval_4(srt_abs_4(1:3),srt_abs_4(1:3))))),'descend');
        
        LeigTemp  = (abs(eigenval_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1)))))^2....
            *(abs(eigenvec_4(1,srt_abs_4(srt_angle_4(1))))^2....
            + abs(eigenvec_4(2,(srt_angle_4(1))))^2....
            + abs(eigenvec_4(4,srt_abs_4(srt_angle_4(1))))^2);
        
        SLeigTemp = (abs(eigenval_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3)))))^2....
            *(abs(eigenvec_4(1,srt_abs_4(srt_angle_4(3))))^2....
            + abs(eigenvec_4(2,srt_abs_4(srt_angle_4(3))))^2....
            + abs(eigenvec_4(4,srt_abs_4(srt_angle_4(3))))^2);
        
        if ((LeigTemp >= SLeigTemp)&&(abs(eigenvec_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1))))...
                >= abs(eigenvec_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3))))) )
            
            vegitation_angle_4(SNR_sample,1) = vegitation_angle_4(SNR_sample,1)....
                + 0.5*angle(eigenval_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3))))/Averaged_samples;
            
            ground_angle_4(SNR_sample,1) = ground_angle_4(SNR_sample,1)....
                + 0.5*angle(eigenval_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1))))/Averaged_samples;
            
            vegitation_abs_4(SNR_sample,1) = vegitation_abs_4(SNR_sample,1)....
                + sqrt(abs(eigenval_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3)))))/Averaged_samples;
                       
            ground_abs_4(SNR_sample,1) = ground_abs_4(SNR_sample,1)....
                + sqrt(abs(eigenval_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1)))))/Averaged_samples;
            
            vegitation_angle_rmse_4(SNR_sample) =  vegitation_angle_rmse_4(SNR_sample)...
                + ((vegitation_offset + angle(eigenval_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3)))))^2)/Averaged_samples;            
            
            ground_angle_rmse_4(SNR_sample) =  ground_angle_rmse_4(SNR_sample)...
                + ((ground_offset - abs(0.5*angle(eigenval_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1))))))^2)/Averaged_samples;            
            
        elseif((LeigTemp >= SLeigTemp)&&(abs(eigenvec_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3))))...
                >= abs(eigenvec_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1))))) )
            
            vegitation_angle_4(SNR_sample,1) = vegitation_angle_4(SNR_sample,1)....
                + 0.5*angle(eigenval_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1))))/Averaged_samples;
            
            ground_angle_4(SNR_sample,1) = ground_angle_4(SNR_sample,1)....
                + 0.5*angle(eigenval_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3))))/Averaged_samples;
            
            vegitation_abs_4(SNR_sample,1) = vegitation_abs_4(SNR_sample,1)....
                + sqrt(abs(eigenval_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1)))))/Averaged_samples;
            
            ground_abs_4(SNR_sample,1) = ground_abs_4(SNR_sample,1)....
                + sqrt(abs(eigenval_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3)))))/Averaged_samples;
            
            vegitation_angle_rmse_4(SNR_sample) =  vegitation_angle_rmse_4(SNR_sample)...
                + ((vegitation_offset + angle(eigenval_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1)))))^2)/Averaged_samples;    
            
            ground_angle_rmse_4(SNR_sample) =  ground_angle_rmse_4(SNR_sample)...
                + ((ground_offset - abs(0.5*angle(eigenval_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3))))))^2)/Averaged_samples;            
            
        elseif (SLeigTemp >= LeigTemp)
            
            vegitation_angle_4(SNR_sample,1) = vegitation_angle_4(SNR_sample,1)....
                + 0.5*angle(eigenval_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1))))/Averaged_samples;
            
            ground_angle_4(SNR_sample,1) = ground_angle_4(SNR_sample,1)....
                + 0.5*angle(eigenval_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3))))/Averaged_samples;
            
            vegitation_abs_4(SNR_sample,1) = vegitation_abs_4(SNR_sample,1)....
                + sqrt(abs(eigenval_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1)))))/Averaged_samples;
            
            ground_abs_4(SNR_sample,1) = ground_abs_4(SNR_sample,1)....
                + sqrt(abs(eigenval_4(srt_abs_4(srt_angle_4(3)),srt_abs_4(srt_angle_4(3)))))/Averaged_samples;
           
            vegitation_angle_rmse_4(SNR_sample) =  vegitation_angle_rmse_4(SNR_sample)...
                + ((vegitation_offset + angle(eigenval_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1)))))^2)/Averaged_samples;        
            
            ground_angle_rmse_4(SNR_sample) =  ground_angle_rmse_4(SNR_sample)...
                + ((ground_offset - abs(0.5*angle(eigenval_4(srt_abs_4(srt_angle_4(1)),srt_abs_4(srt_angle_4(1))))))^2)/Averaged_samples;            
        end
    end
end
%% Plotting Results

figure(1);
title('4rth Order ESPRIT Interferometric Phases');
xlabel('SNR (dB)');ylabel('Interferometric Phase (degrees)');
hold on;
plot(SNR,(-ground_angle_4)*180/pi,'bx');
plot(SNR,(-vegitation_angle_4)*180/pi,'gx');
plot(SNR,G_O*ones(1,SNR_samples),'b');
plot(SNR,V_O*ones(1,SNR_samples),'g');
axis([-15,15,G_O-5,V_O+5])
legend('4rth Order Ground Estimate','4rth Order Vegetation Estimate',...
    '4rth Order Ground Actual','4rth Order Vegetation Actual','Location','east')
hold off

figure(3);
title('4rth Order ESPRIT Magnitude Plot');
xlabel('SNR (dB)');ylabel('Magnitude');
hold on;
plot(SNR,ground_abs_4,'bx');
plot(SNR,vegitation_abs_4,'gx');
% axis([-15,15,0,2])
legend('4rth Order Ground','4rth Order Vegetation','Location','northwest')
hold off

figure(5)
title('4rth Order ESPRIT RMS Error Vs. SNR');
xlabel('SNR (dB)');ylabel('RMS Error (degrees)');
hold on;
plot(SNR,ground_angle_rmse_4,'bx')
plot(SNR,vegitation_angle_rmse_4,'gx')
legend('4rth Order Ground','4rth Order Vegetation','Location','northwest')
hold off

figure(7)
title('4rth Order ESPRIT Magnitude Vs.RMS Error');
ylabel('RMS Error (degrees)');xlabel('Magnitude');
hold on;
plot(ground_abs_4,ground_angle_rmse_4,'bx')
plot(vegitation_abs_4,vegitation_angle_rmse_4,'gx')
legend('4rth Order Ground','4rth Order Vegetaion','Location','northwest');
hold off