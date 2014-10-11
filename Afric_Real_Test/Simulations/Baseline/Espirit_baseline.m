%% Espirit_baseline
clear;clc;
%% Initializations
g_weight = 1;
v_weight = 0.1;
n_weight = 0;
g_angle = pi/6;
v_angle = pi/4;
array_size = 6;
Averaged_samples = 1000;

X = zeros(1,array_size);

phase = zeros(1,Averaged_samples);

%% Two Element Array size
% A on Case 1

for Averaged_sample = 1:Averaged_samples;
    
    phi_one = 1i*2*pi*rand;
    phi_two = 1i*2*pi*rand;
    
    for pos = 1:array_size;
        
        distance = 2*pi*pos/2;  % element separation of 1/2 wavelength
        
        X(1,pos) = g_weight*exp(phi_one)*exp(-1i*distance*sin(g_angle))... %Ground
            + v_weight*exp(phi_two)*exp(-1i*distance*sin(v_angle))...  %Vegetaion
            + n_weight*sqrt(-2.*log(1-rand))*exp(1i*2*pi*rand);  %Noise
    end
    
    Y1 = [X(1,1); X(1,2)]; Y2 = [X(1,5);X(1,6)]; %Baseline = 4d
    
    % Y1 = [X(1,2),X(1,3)]; Y2 = [X(1,4),X(1,5)]; %Baseline = 2d
    %    Y1 = [X(1,3),X(1,4)]; Y2 = [X(1,4),X(1,5)]; %Baseline = d
    %     Y1 = [X(1,2),X(1,3),X(1,4)]; Y2 = [X(1,3),X(1,4),X(1,5)]; %Baseline = d
    %     Y1 = [X(1,1),X(1,2),X(1,3),X(1,4)]; Y2 = [X(1,2),X(1,3),X(1,4),X(1,5)]; %Baseline = d
    
    R1 = bsxfun(@times, Y1,Y1'); % bsxfun is a Faster outer Product than Y1*Y1'
    R2 =  bsxfun(@times, Y1,Y2');
    A =  bsxfun(@times,pinv(R1),R2);
    
    [u,uv] = eig(A);
    [~,kk]=sort(angle(diag(uv)),'ascend');
    
    phase(1,Averaged_sample) = phase(Averaged_sample) + angle(uv(kk(1),kk(1)));
  
    
end

%% Ploting Results

figure(1)
plot(phase,'b+');

