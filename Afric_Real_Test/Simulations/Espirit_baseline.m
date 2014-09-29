%% Afric_Sim_test_Cumulant
clear;clc;
%% Initializations
mag_one = 1;
mag_two = 1/2;
phase_one = pi/6;
phase_two = pi/4;
array_size = 6;
Averaged_samples = 1000;
Noise_weight = 0.1;
X = zeros(1,array_size);
phi1=zeros(Averaged_samples,1);
mag1=zeros(Averaged_samples,1);


%% Two Element Array size
% A on Case 1
for Averaged_sample = 1:Averaged_samples;
    %     Noise = Averaged_sample*Noise_weight;
    for pos = 1:array_size;
        X(1,pos) = mag_one*exp(1i*2*pi*rand)*exp(-1i*sin(phase_one))...
            + mag_two*exp(1i*2*pi*rand)*exp(-1i*sin(phase_two))...
            + sqrt(-2*log(1-rand)).*exp(1i*2*pi*rand);
    end
    
     Y1 = [X(1,1),X(1,2)]; Y2 = [X(1,5),X(1,6)]; %Baseline = 4d
%     Y1 = [X(1,2),X(1,3)]; Y2 = [X(1,4),X(1,5)]; %Baseline = 2d
%     Y1 = [X(1,3),X(1,4)]; Y2 = [X(1,4),X(1,5)]; %Baseline = d
%     Y1 = [X(1,2),X(1,3),X(1,4)]; Y2 = [X(1,3),X(1,4),X(1,5)]; %Baseline = d
%     Y1 = [X(1,1),X(1,2),X(1,3),X(1,4)]; Y2 = [X(1,2),X(1,3),X(1,4),X(1,5)]; %Baseline = d
    
    R1 = Y1*Y1';
    R2 = Y1*Y2';
    
    A = R1*R2';
   
    phi1(Averaged_sample)=phi1(Averaged_sample) + angle(A);
    mag1(Averaged_sample)=mag1(Averaged_sample) + abs(A);
    
end
%% Ploting Results

figure(1)
plot(phi1*180/pi,'b+');

figure(2)
plot(mag1,'b+');
