%%Espirit Algorithm Fourth order
clear;clc;
%% Initializations
fid(1)=fopen('ref_hh.dat','r');fid(2)=fopen('ref_vv.dat','r');
fid(3)=fopen('ref_xx.dat','r');fid(4)=fopen('offset_hh.dat','r');
fid(5)=fopen('offset_vv.dat','r');fid(6)=fopen('offset_xx.dat','r');

h1=fread(fid(1),[730 1662],'single');h2=fread(fid(4),[730 1662],'single');
v1=fread(fid(2),[730 1662],'single');v2=fread(fid(5),[730 1662],'single');
x1=fread(fid(3),[730 1662],'single');x2=fread(fid(6),[730 1662],'single');

for close = 1:6;
fclose(fid(close));
end

hh.ref=h1(1:2:729,:) + 1i*h1(2:2:730,:); hh.off=h2(1:2:729,:) + 1i*h2(2:2:730,:);
vv.ref=v1(1:2:729,:) + 1i*v1(2:2:730,:); vv.off=v2(1:2:729,:) + 1i*v2(2:2:730,:);
xx.ref=x1(1:2:729,:) + 1i*x1(2:2:730,:); xx.off=x2(1:2:729,:) + 1i*x2(2:2:730,:);

hhref = hh.ref;hhoff = hh.off;
vvref = vv.ref;vvoff = vv.off;
xxref = xx.ref;xxoff = xx.off;

[ylength,xlength] = size(hhoff);
r = 7; c = 7;
L = r+c+1;
Lreshape = (r+c+1)^2;
eye_w = 0.3737;
g = zeros(ylength-r,xlength-c);
v = zeros(ylength-r,xlength-c);

forth_order = 6;
sort_ground_4 = zeros(forth_order,1);
sort_vegetation_4 = zeros(forth_order,1);
%% Cumulant Martix Calculations
%%%%%%%%%%%%%%%%%%%%%Pull rand line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = r+1:ylength-r;
    for col = c+1:xlength-c;
        
        Var(1:L,1:L) = hhref(row-r:row+r,col-c:col+c);
        sh1(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var(1:L,1:L) = hhoff(row-r:row+r,col-c:col+c);
        sh2(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var(1:L,1:L) = vvref(row-r:row+r,col-c:col+c);
        sv1(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var(1:L,1:L) = vvoff(row-r:row+r,col-c:col+c);
        sv2(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var(1:L,1:L) = xxref(row-r:row+r,col-c:col+c);
        sx1(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var(1:L,1:L) = xxoff(row-r:row+r,col-c:col+c);
        sx2(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        S1 = [sh1.*sh1;sv1.*sv1;sx1.*sx1;sh1.*sv1;sh1.*sx1;sv1.*sx1;];
        S2 = [sh2.*sh2;sv2.*sv2;sx2.*sx2;sh2.*sv2;sh2.*sx2;sv2.*sx2;];
        
        R1_4 = S1*S1';  
        
        R2_4 = S1*S2';        
        
        [eigenvec_4,eigenval_4] = eig(pinv(R1_4 + eye_w*eye(6))*R2_4);
        
        for i = 1:forth_order
            sort_ground_4(i) = (abs(eigenval_4(i,i)))^2*(abs(eigenvec_4(1,i))^2 + abs(eigenvec_4(2,i))^2 + abs(eigenvec_4(4,i))^2);
            sort_vegetation_4(i) = (abs(eigenval_4(i,i)))^2*(abs(eigenvec_4(3,i))^2 + abs(eigenvec_4(5,i))^2 + abs(eigenvec_4(6,i))^2);
        end       
        
        [~,srt_g_4] = sort(sort_ground_4,'descend');
        [~,srt_v_4] = sort(sort_vegetation_4,'descend');
        
        g(row,col) = eigenval_4(srt_g_4(1),srt_g_4(1));
        v(row,col) = eigenval_4(srt_v_4(1),srt_v_4(1));
        
    end
end
%% Plotting gv Results
figure(1); imagesc(0.5*angle(g)); title('4th Ord angle(g)');
figure(2); imagesc(0.5*angle(v)); title('4th Ord angle(v)');
figure(3); imagesc(abs(g)); title('4th Ord abs(g)');
figure(4); imagesc(abs(v)); title('4th Ord abs(v)');