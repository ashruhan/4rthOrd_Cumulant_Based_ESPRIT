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

one_4 = zeros(ylength-r,xlength-c);
two_4 = zeros(ylength-r,xlength-c);
three_4 = zeros(ylength-r,xlength-c);
four_4 = zeros(ylength-r,xlength-c);
five_4 = zeros(ylength-r,xlength-c);
six_4 = zeros(ylength-r,xlength-c);

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

        [~,srt_4] = sort(abs(eigenval_4),'descend');
        
        one_4(row,col) = eigenval_4(srt_4(1),srt_4(1));
        two_4(row,col) = eigenval_4(srt_4(2),srt_4(2));
        three_4(row,col) = eigenval_4(srt_4(3),srt_4(3));
        four_4(row,col) = eigenval_4(srt_4(4),srt_4(4));
        five_4(row,col) = eigenval_4(srt_4(5),srt_4(5));
        six_4(row,col) = eigenval_4(srt_4(6),srt_4(6));
        
    end
end
%% Plotting gv Results
