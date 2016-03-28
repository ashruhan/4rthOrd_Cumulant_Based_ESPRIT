%% Afric_Sim_test_ESPRIT
clc;clear;
%% Initializations
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor
Averaged_samples = 100;
Window_optimal = 81;    %size of window

Pol_ground = [1;1;0]/sqrt(2);
Pol_Cum_ground = [1;1;0;1;0;0]/sqrt(3); %ground
Pol_vegitation = [1;1;1]/sqrt(3);
Pol_Cum_vegitation = [1;1;1;1;1;1]/sqrt(6); %vegitation

SNR_Distribution = 31;
SNR = linspace(-15,15,SNR_Distribution);

G_O = 20;
ground_offset = G_O*pi/180; % ground interferomitry offset

VO_DistLength = 80;
V_O_Dist = linspace(10,89,VO_DistLength);
vegetation_offset = V_O_Dist.*pi/180;    % veg interferomitry offset

eye_distlength = 50;
eye_2 = linspace(-1,1,eye_distlength);
eye_4 = linspace(-1,1,eye_distlength);

ground_angle_4 = zeros(SNR_Distribution,VO_DistLength,eye_distlength);
ground_angle_2 = zeros(SNR_Distribution,VO_DistLength,eye_distlength);

vegitation_angle_4 = zeros(SNR_Distribution,VO_DistLength,eye_distlength);
vegitation_angle_2 = zeros(SNR_Distribution,VO_DistLength,eye_distlength);

ground_abs_4 = zeros(SNR_Distribution,VO_DistLength,eye_distlength);
ground_abs_2 = zeros(SNR_Distribution,VO_DistLength,eye_distlength);

vegitation_abs_4 = zeros(SNR_Distribution,VO_DistLength,eye_distlength);
vegitation_abs_2 = zeros(SNR_Distribution,VO_DistLength,eye_distlength);

%% Matrix Construction
for SNR_sample = 1:SNR_Distribution
    Noise = (10^(-SNR(SNR_sample)/20))/sqrt(3);
    for V_O_index = 1:VO_DistLength
        for eye_index = 1:eye_distlength
            for unusedvariable = 1:Averaged_samples
                
                %% Matrix Construction
                
                g =  Pol_ground*exp(1i*2*pi*rand(1,Window_optimal));
                v =  Pol_vegitation*exp(1i*2*pi*rand(1,Window_optimal));
                
                s1 = alpha*g + beta*v;
                s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*vegetation_offset(V_O_index))*v;
                
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
                
                [eigenvec_2,eigenval_2] =eig((pinv(R1_2 + eye_2(eye_index)*eye(3)))...
                    *R2_2,'nobalance');
                
                polfilter_2 = abs(Pol_ground'*eigenvec_2);
                [~,srt_2] = sort(polfilter_2,'descend');
                ground_angle_2(SNR_sample,V_O_index,eye_index) = ground_angle_2(SNR_sample,V_O_index,eye_index)...
                    + ((ground_offset - abs(angle(eigenval_2(srt_2(1),srt_2(1)))))^2)/Averaged_samples;
                
                ground_abs_2(SNR_sample,V_O_index,eye_index) = ground_abs_2(SNR_sample,V_O_index,eye_index)...
                    + abs(eigenval_2(srt_2(1),srt_2(1)))/Averaged_samples;
                
                polfilter_2 = abs(Pol_vegitation'*eigenvec_2);
                [~,srt_2] = sort(polfilter_2,'descend');
                vegitation_angle_2(SNR_sample,V_O_index,eye_index) = vegitation_angle_2(SNR_sample,V_O_index,eye_index)...
                    + ((vegetation_offset(V_O_index) - abs(angle(eigenval_2(srt_2(1),srt_2(1)))))^2)/Averaged_samples;
                
                vegitation_abs_2(SNR_sample,V_O_index,eye_index) = vegitation_abs_2(SNR_sample,V_O_index,eye_index)...
                    + abs(eigenval_2(srt_2(1),srt_2(1)))/Averaged_samples;
                
                %% Fourth Order Statistics
                [ Cumulant_11, Cumulant_12] = Cumulant( s1_Noise ,s2_Noise,Window_optimal );
                
                [~,eigenvalCov_4] = eig(Cumulant_11,'nobalance');
                
                [eigenvec_4,eigenval_4] = eig((pinv(Cumulant_11 + eye_4(eye_index)*eye(6)))...
                    *Cumulant_12,'nobalance');
                
                polfilter_4 = abs(Pol_Cum_ground'*eigenvec_4);
                [~,srt_4] = sort(polfilter_4,'descend');
                ground_angle_4(SNR_sample,V_O_index,eye_index) = ground_angle_4(SNR_sample,V_O_index,eye_index)...
                    + ((ground_offset - abs(0.5*angle(eigenval_4(srt_4(1),srt_4(1)))))^2)/Averaged_samples;
                
                ground_abs_4(SNR_sample,V_O_index,eye_index) = ground_abs_4(SNR_sample,V_O_index,eye_index)...
                    + abs(eigenval_4(srt_4(1),srt_4(1)))/Averaged_samples;
                
                polfilter_4 = abs(Pol_Cum_vegitation'*eigenvec_4);
                [~,srt_4] = sort(polfilter_4,'descend');
                vegitation_angle_4(SNR_sample,V_O_index,eye_index) = vegitation_angle_4(SNR_sample,V_O_index,eye_index)...
                    + ((vegetation_offset(V_O_index) - abs(0.5*angle(eigenval_4(srt_4(1),srt_4(1)))))^2)/Averaged_samples;
                
                vegitation_abs_4(SNR_sample,V_O_index,eye_index) = vegitation_abs_4(SNR_sample,V_O_index,eye_index)...
                    + abs(eigenval_4(srt_4(1),srt_4(1)))/Averaged_samples;
            end
        end
    end
end

ground_angle_4 = sqrt(ground_angle_4)*180/pi;
vegitation_angle_4 = sqrt(vegitation_angle_4)*180/pi;
ground_angle_2 = sqrt(ground_angle_2)*180/pi;
vegitation_angle_2 = sqrt(vegitation_angle_2)*180/pi;
%% Permutate
ground_angle_4 = permute(ground_angle_4,[2 3 1]);
vegitation_angle_4 = permute(vegitation_angle_4,[2 3 1]);
ground_angle_2 = permute(ground_angle_2,[2 3 1]);
vegitation_angle_2 = permute(vegitation_angle_2,[2 3 1]);

ground_abs_4 = permute(ground_abs_4,[2 3 1]);
vegitation_abs_4 = permute(vegitation_abs_4,[2 3 1]);
ground_abs_2 = permute(ground_abs_2,[2 3 1]);
vegitation_abs_2 = permute(vegitation_abs_2,[2 3 1]);

%% Plotting Results
SNRPLOT = SNR_Distribution-10;
figure(1);

meshc(ground_angle_4(:,:,SNRPLOT));title('4th order ground RMS Error');
xlabel('vegetation sweep 1-89 degrees');ylabel('ground sweep 1-89 degrees');
figure(2);
meshc(vegitation_angle_4(:,:,SNRPLOT));title('4th order vegetation RMS Error');
xlabel('vegetation sweep 1-89 degrees');ylabel('ground sweep 1-89 degrees');
figure(3)
meshc(ground_angle_2(:,:,SNRPLOT));title('2nd order ground RMS Error');
xlabel('vegetation sweep 1-89 degrees');ylabel('ground sweep 1-89 degrees');
figure(4)
meshc(vegitation_angle_2(:,:,SNRPLOT));title('2nd order vegetation RMS Error');
xlabel('vegetation sweep 1-89 degrees');ylabel('ground sweep 1-89 degrees');
figure(5)
meshc(ground_abs_4(:,:,SNRPLOT));title('4th order ground coherance');
xlabel('vegetation sweep 1-89 degrees');ylabel('ground sweep 1-89 degrees');
figure(6)
meshc(vegitation_abs_4(:,:,SNRPLOT));title('4th order vegetation coherance');
xlabel('vegetation sweep 1-89 degrees');ylabel('ground sweep 1-89 degrees');
figure(7)
meshc(ground_abs_2(:,:,SNRPLOT));title('2nd order ground coherance');
xlabel('vegetation sweep 1-89 degrees');ylabel('ground sweep 1-89 degrees');
figure(8)
meshc(vegitation_abs_2(:,:,SNRPLOT));title('2nd order vegetation coherance');
xlabel('vegetation sweep 1-89 degrees');ylabel('ground sweep 1-89 degrees');