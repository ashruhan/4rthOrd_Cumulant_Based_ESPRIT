%% Espirit_baseline
clear;clc;
%% Initializations
ground = 1;
vegitation = 1/2;
noiseweight = 0.1;
gphse = pi/6;
vphase = pi/4;
array_size = 6;
Averaged_samples = 1000;

X = zeros(1,array_size);
phase1=zeros(Averaged_samples,1);
phase2=zeros(Averaged_samples,1);


%% Two Element Array size
% A on Case 1
for Averaged_sample = 1:Averaged_samples;
    phi_one = rand;
    phi_two = rand;
    
    
    for pos = 1:array_size;
        
        distance = pos/2;
        X(1,pos) = ground*exp(1i*2*pi*phi_one)*exp(-1i*distance*sin(gphse))... %Ground
            + vegitation*exp(1i*2*pi*phi_two)*exp(-1i*distance*sin(vphase))...  %Vegetaion
            + noiseweight*sqrt(-2*log(1-rand)).*exp(1i*2*pi*rand);  %Noise
    end
    
    %     Y1 = [X(1,1),X(1,2)]; Y2 = [X(1,5),X(1,6)]; %Baseline = 4d
   % Y1 = [X(1,2),X(1,3)]; Y2 = [X(1,4),X(1,5)]; %Baseline = 2d
    %    Y1 = [X(1,3),X(1,4)]; Y2 = [X(1,4),X(1,5)]; %Baseline = d
    %     Y1 = [X(1,2),X(1,3),X(1,4)]; Y2 = [X(1,3),X(1,4),X(1,5)]; %Baseline = d
    %     Y1 = [X(1,1),X(1,2),X(1,3),X(1,4)]; Y2 = [X(1,2),X(1,3),X(1,4),X(1,5)]; %Baseline = d
    
    R1 = bsxfun(@times, Y1,Y1');
    R2 =  bsxfun(@times, Y1,Y2');
    
    A =  bsxfun(@times,R1,R2');
    
    [u,uv] = eig(A);
    [~,kk]=sort(angle(diag(uv)),'ascend');
    
    phase1(Averaged_sample)=phase1(Averaged_sample) + angle(uv(kk(1),kk(1)));
    phase2(Averaged_sample)=phase2(Averaged_sample) + angle(uv(kk(2),kk(2)));
    
end
%% Ploting Results

figure(1)
plot(phase1,'b+');

figure(2)
plot(phase2,'b+');
