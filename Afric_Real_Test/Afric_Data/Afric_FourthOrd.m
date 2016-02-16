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

g = zeros(ylength-r,xlength-c);
v = zeros(ylength-r,xlength-c);

n1 = zeros(ylength-r,xlength-c);
n2 = zeros(ylength-r,xlength-c);
n3 = zeros(ylength-r,xlength-c);
n4 = zeros(ylength-r,xlength-c);

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
        
        R1 = S1*S1';  R2 = S1*S2';
   
        [eigenvec,eigenval] = eig(pinv(R1)*R2);

        [~,srt]=sort(angle(diag(eigenval)),'descend');
        
        Leig_copol = abs(eigenvec(1,srt(1))) + abs(eigenvec(2,srt(1)));
        SLeig_copol = abs(eigenvec(1,srt(2))) + abs(eigenvec(2,srt(2)));
        
        n1(row,col) = eigenval(srt(3),srt(3));
        n2(row,col) = eigenval(srt(4),srt(4));
        n3(row,col) = eigenval(srt(5),srt(5));
        n4(row,col) = eigenval(srt(6),srt(6));
        
        if (Leig_copol >= SLeig_copol)
            
           g(row,col) = eigenval(srt(1),srt(1));
           v(row,col) = eigenval(srt(2),srt(2));
           
        else    
            
           g(row,col) = eigenval(srt(2),srt(2));
           v(row,col) = eigenval(srt(1),srt(1));
           
        end
        
    end
end
%% Plotting gv Results
figure(1); imagesc(0.5*angle(g)); title('4th Ord angle(g)');
figure(2); imagesc(0.5*angle(v)); title('4th Ord angle(v)');
figure(3); imagesc(abs(g)); title('4th Ord abs(g)');
figure(4); imagesc(abs(v)); title('4th Ord abs(v)');

%% Plotting Noise Results
figure(5); imagesc(0.5*angle(n1)); title('4th Ord angle(n1)');
figure(6); imagesc(0.5*angle(n2)); title('4th Ord angle(n2)');
figure(7); imagesc(0.5*angle(n3)); title('4th Ord angle(n3)');
figure(8); imagesc(0.5*angle(n4)); title('4th Ord angle(n4)');

figure(9); imagesc(abs(n1)); title('4th Ord abs(n1)');
figure(10); imagesc(abs(n2)); title('4th Ord abs(n2)');
figure(11); imagesc(abs(n3)); title('4th Ord abs(n3)');
figure(12); imagesc(abs(n4)); title('4th Ord abs(n4)');