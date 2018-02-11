%%%%%%%%%%%%%%%%%%%%%%%%%%%%DECLARATIONS%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
pg1 = [1;0;1]/sqrt(2); %ground 
pv1 = [1;1;1]/sqrt(3); %veg
ground_offset = pi/2; % ground interferomitry offset 
veg_offset = pi/6;    % veg interferomitry offset

Window = 100;    %size of window
average = 100;  %Samples
  

for Noise = 100:-1:1;

   NN(Noise)=1./(0.01*Noise)^2;
   RMS_veg = 0;    %used for averaging
   RMS_ground = 0;  %used for averaging  
    
for Ave = 1:average;

g =  pg1*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
v =  pv1*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));

s1 = g + v;
s2 = exp(1i*ground_offset)*g + exp(1i*veg_offset)*v;

s1_Noise = s1 + .01*Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
s2_Noise = s2 + .01*Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
 
 %%%%%%%%%%%%%%%%%%Super Resolution%%%%%%%%%%%%%%%%%%%%%%
 R1 = s1_Noise*s1_Noise';
 R2 = s2_Noise*s1_Noise';
 [u ,c] = eig(R2*R2');
 R1_bar = u*sqrt(c)*u';
 Q = pinv(R1_bar)*R2;
 R2_bar = R1_bar*Q*Q;
  %%%%%%%%%%%%%%%%%%Super Resolution%%%%%%%%%%%%%%%%%%%%%%

 [U ,Uv] = eig(R2_bar*pinv(R1));
 %%%%%%%%%%%%%%%%%%%%%%%%%Sorting out Ground from Vege%%%%%%%%%%%%%%%%%%%%
 
 %%%SORTING ROUTING WILL ONLY WORK WHEN THE INTERFEROMIC PHASE DIFFERANCE OF
 %%%GROUND IS LARGER THAN THE PHASE DIFFERANCE OF THE VEGITATION
[UU,kk]=sort(abs(angle(diag(Uv))-pi/2),'ascend');
GG = Uv(kk(1),kk(1));

[UU,kk]=sort(abs(angle(diag(Uv))-pi/6),'ascend');
VV = Uv(kk(1),kk(1));


 %%%%%%%%%%%%%%%%%%%%%%%End Sorting out Ground from Vege%%%%%%%%%%%%%%%%%%%
 
RMS_ground = RMS_ground + (angle(GG)-ground_offset)^2; %%SUMMING OVER "AVERAGE" AMOUNT OF POINTS
RMS_veg = RMS_veg +(angle(VV)-veg_offset)^2; %%SUMMING OVER "AVERAGE" AMOUNT OF POINTS

end

RMS_ground = sqrt((1/(2*average))*RMS_ground);%%FINAL AVERAGING
RMS_veg = sqrt((1/(2*average))*RMS_veg);      %%FINAL AVERAGING

RMS_G(Noise) = RMS_ground; %%SAVING VALUES FOR A SET POINT OF NOISE THEN PLOTTING AFTER
RMS_V(Noise) = RMS_veg;    %%SAVING VALUES FOR A SET POINT OF NOISE THEN PLOTTING AFTER

end

 figure 
 hold on
 plot(10*log10(NN),180*RMS_G/pi)
 plot(10*log10(NN),180*RMS_V/pi,'r');
 xlabel('S/N Ratio, dB'); ylabel('Error In Degrees');legend('Ground','Vegetation','NorthEast')
 hold off



 
 