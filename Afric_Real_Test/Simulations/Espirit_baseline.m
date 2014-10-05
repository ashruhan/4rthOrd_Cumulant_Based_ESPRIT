%% Espirit_baseline
clear;clc;
%% Initializations
g_weight = 1;
v_weight = 0.5;
n_weight = 0;
g_phse = pi/6;
v_phase = pi/4;
array_size = 6;
Averaged_samples = 10;

X = zeros(1,array_size);
phase1 = zeros(1,Averaged_samples);
phase2 = zeros(1,Averaged_samples);

comp_1 = zeros(1,Averaged_samples);
comp_2 = zeros(1,Averaged_samples);
%% Two Element Array size
% A on Case 1
for Averaged_sample = 1:Averaged_samples;
    phi_one = rand();
    phi_two = rand();
    for pos = 1:array_size;
        
        distance = pos/2;
        X(1,pos) = g_weight*exp(1i*2*pi*phi_one)*exp(-1i*distance*sin(g_phse))... %Ground
            + v_weight*exp(1i*2*pi*phi_two)*exp(-1i*distance*sin(v_phase))...  %Vegetaion
            + n_weight*sqrt(-2*log(1-rand))*exp(1i*2*pi*rand);  %Noise
    end
    
    Y1 = [X(1,1); X(1,2)]; Y2 = [X(1,5);X(1,6)]; %Baseline = 4d
    % Y1 = [X(1,2),X(1,3)]; Y2 = [X(1,4),X(1,5)]; %Baseline = 2d
    %    Y1 = [X(1,3),X(1,4)]; Y2 = [X(1,4),X(1,5)]; %Baseline = d
    %     Y1 = [X(1,2),X(1,3),X(1,4)]; Y2 = [X(1,3),X(1,4),X(1,5)]; %Baseline = d
    %     Y1 = [X(1,1),X(1,2),X(1,3),X(1,4)]; Y2 = [X(1,2),X(1,3),X(1,4),X(1,5)]; %Baseline = d
    
    R1 = bsxfun(@times, Y1,Y1');
    R2 =  bsxfun(@times, Y1,Y2');
    
    A =  bsxfun(@times,R1,R2');
    
    [u,uv] = eig(A);
    [~,kk]=sort(angle(diag(uv)),'ascend');
    
    phase1(1,Averaged_sample) = phase1(Averaged_sample) + angle(uv(kk(1),kk(1)));
    phase2(1,Averaged_sample) = phase2(Averaged_sample) + angle(uv(kk(2),kk(2)));
    comp_1(1,Averaged_sample) = comp_1(Averaged_sample) + phi_one;
    comp_2(1,Averaged_sample) = comp_2(Averaged_sample) + phi_two;
end

%% Ploting Results

figure(1)
plot(phase1,'b+');
hold on;
plot(comp_1)
hold off;

figure(2)
plot(phase2,'b+');hold on;
plot(comp_2);
hold off;