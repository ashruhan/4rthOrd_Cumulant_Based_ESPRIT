alpha = 1; %ground weighting factor
beta = 1.0;   %veg weighting factor
pg1 = [1;1;0]; %ground 
pv1 = [1;1;1]; %veg
ground_offset = -pi/3; % ground interferomitry offset 
veg_offset = pi/6;    % veg interferomitry offset
x = .001:100:1; %% used for graphing
Window = 100;    %size of window
N0=100;

phi1=zeros(10,1);
phi2=zeros(10,1);
mag1=phi1;
mag2=phi2;

for in=10:10
    
Noise=in*1;

for it=1:N0
    
g =  pg1*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));
v =  pv1*(sqrt(-2*log(1-rand(1,Window))).*exp(1i*2*pi*rand(1,Window)));

s1 = alpha*g + beta*v;
s2 = alpha*exp(1i*ground_offset)*g + beta*exp(1i*veg_offset)*v;

s1_Noise = s1 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
s2_Noise = s2 + Noise*sqrt(-2*log(1-rand(3,Window))).*exp(1i*2*pi*rand(3,Window));
 
 R1 = s1_Noise*s1_Noise';
 R2 = s1_Noise*s2_Noise';
 A=R2*pinv(R1);
 
 [u,c] = eig(A^0.5);
 
 s1=abs(pg1'*u);
 [ss,kk]=sort(s1,'descend');
 phi1(in)=phi1(in) + 2*angle(c(kk(1),kk(1)))/N0;
 mag1(in)=mag1(in) + abs(c(kk(1),kk(1)))/N0;

 s2=abs(pv1'*u);
 [ss,kk]=sort(s2,'descend');
 phi2(in)=phi2(in) + 2*angle(c(kk(1),kk(1)))/N0;
 mag2(in)=mag2(in) + abs(c(kk(1),kk(1)))/N0;

end
end
zplane(eig(A));
figure
hold
plot(phi1*180/pi,'b+');
plot(phi2*180/pi,'ro');
% 
% figure
% hold
% plot(mag1,'b+');
% plot(mag2,'ro');
