%%%%%%%%%%%%%%%%%%%%%%%%%%%%DECLARATIONS%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
alpha = 1; %ground weighting factor
beta = 1;   %veg weighting factor
pg1 = [1;0;1]/sqrt(2); %ground
pv1 = [1;1;1]/sqrt(3); %veg
ground_offset = pi/2; % ground interferomitry offset
veg_offset = pi/6;    % veg interferomitry offset

Window = 25;    %size of window
average = 100;  %Samples

Eig_R1_s = zeros(1,100);
Eig_R1_m = zeros(1,100);
Eig_R2_s = zeros(1,100);
Eig_R2_m = zeros(1,100);
N_R1 = zeros(1,100);
N_R2 = zeros(1,100);
%%%%%%%%%%%%%%%%%%%%%%%%End DECLARATIONS%%%%%%%%%%%%%%%%%%%%%%
for Noise = 100:-1:1;
    
    NN(Noise)=1./(0.01*Noise)^2;

    g =  pg1*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
    v =  pv1*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
    
    s1 = alpha*g + beta*v;
    s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*veg_offset)*v;
    
    s1_Noise = s1 + .01*Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
    s2_Noise = s2 + .01*Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
    
    %%%%%%%%%%%%%%%%%%Super Resolution%%%%%%%%%%%%%%%%%%%%%%
    R1 = s1_Noise*s1_Noise';
    R2 = s2_Noise*s1_Noise';
    [u ,c] = eig(R2'*R2);
    
    H = u*sqrt(c)*u';
    
    C = R2*H';
    Z = H*R1';
    %%%%%%%%%%%%%%%%%%Super Resolution%%%%%%%%%%%%%%%%%%%%%%
    
    [Uc ,Uvc] = eig(C);
    [Uz ,Uvz] = eig(Z);
    
    [Sort_c,kc]=sort(diag(Uc),'ascend');
    [Sort_z,kz]=sort(diag(Uz),'ascend');
    
    Eig_R2_s(Noise) = 1- Sort_c(1,1);
    Eig_R2_m(Noise) = 1- Sort_c(2,1);
    Eig_R1_s(Noise) = 1- Sort_z(1,1);
    Eig_R1_m(Noise) = 1- Sort_z(2,1);
    N_R2(Noise) = 1- Sort_c(3,1);
    N_R1(Noise) = 1- Sort_z(3,1);

end

 figure(1)
subplot(4,1,1)
plot(10*log10(NN),(180/pi)*angle(Eig_R1_s));xlabel('S/N Ratio, dB'); ylabel('Error In Degrees');
subplot(4,1,2)
plot(10*log10(NN),(180/pi)*angle(Eig_R1_m));xlabel('S/N Ratio, dB'); ylabel('Error In Degrees');

subplot(4,1,3)
plot(10*log10(NN),abs(Eig_R1_s));xlabel('S/N Ratio, dB'); ylabel('Error In Magnitude');
subplot(4,1,4)
plot(10*log10(NN),abs(Eig_R1_m));xlabel('S/N Ratio, dB'); ylabel('Error In Magnitude');
hold on;

figure(2)
subplot(4,1,1)
plot(10*log10(NN),(180/pi)*angle(Eig_R2_s));xlabel('S/N Ratio, dB'); ylabel('Error In Degrees');
subplot(4,1,2)
plot(10*log10(NN),(180/pi)*angle(Eig_R2_m));xlabel('S/N Ratio, dB'); ylabel('Error In Degrees');

subplot(4,1,3)
plot(10*log10(NN),abs(Eig_R2_s));xlabel('S/N Ratio, dB'); ylabel('Error In Magnitude');
subplot(4,1,4)
plot(10*log10(NN),abs(Eig_R2_m));xlabel('S/N Ratio, dB'); ylabel('Error In Magnitude');
hold off;
