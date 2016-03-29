%% Espirit Algorithm Second order
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

r = 7; c = 7;
L = r+c+1;
Lreshape = (r+c+1)^2;
[ylength,xlength] = size(hhoff);
ground_2 = zeros(ylength-r,xlength-c);
vegetation_2 = zeros(ylength-r,xlength-c);
%% Espirit Martix Calculations
%%%%%%%%%%%%%%%%%%%%%Pull rand line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = fliplr(r+1:ylength-r);
    for col = c+1:xlength-c;
        
        Var = hhref(row-r:row+r,col-c:col+c);
        S1(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var = hhoff(row-r:row+r,col-c:col+c);
        S2(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var = vvref(row-r:row+r,col-c:col+c);
        S1(2,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var = vvoff(row-r:row+r,col-c:col+c);
        S2(2,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var = xxref(row-r:row+r,col-c:col+c);
        S1(3,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var = xxoff(row-r:row+r,col-c:col+c);
        S2(3,1:Lreshape) = reshape(Var,1,Lreshape);
        
        R1 = S1*S1';  R2 = S1*S2';
        
        [eigenvec,eigenval] = eig(pinv(R1)*R2);
        
        [~,srt_2]=sort(abs(diag(eigenval)),'descend');
        
        Leig_copol = abs(eigenvec(1,srt_2(1)))^2 + abs(eigenvec(2,srt_2(1)))^2;
        SLeig_copol = abs(eigenvec(1,srt_2(2)))^2 + abs(eigenvec(2,srt_2(2)))^2;
                
        if (Leig_copol >= SLeig_copol)
            
           ground_2(row,col) = eigenval(srt_2(1),srt_2(1));
           vegetation_2(row,col) = eigenval(srt_2(2),srt_2(2));
           
        else    
            
           ground_2(row,col) = eigenval(srt_2(2),srt_2(2));
           vegetation_2(row,col) = eigenval(srt_2(1),srt_2(1));
           
        end
        
    end
end

%% Plotting Results
figure(1); imagesc(angle(ground_2)); title('2nd Order Ground Interferometric Phase');
figure(2); imagesc(angle(vegetation_2)); title('2nd Order Vegetation Interferometric Phase');
figure(3); imagesc(angle(vegetation_2) - angle(ground_2)); title('4th Order V-G');
figure(4); imagesc(abs(ground_2)); title('2nd Order Ground Magnitude');
figure(5); imagesc(abs(vegetation_2)); title('2nd Order Vegetation Magnitude');