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

halfrow = 7; halfcol = 7;
base = halfrow+halfcol+1;
Window_optimal = (base)^2;

ground_4 = zeros(ylength-halfrow,xlength-halfcol);
vegetation_4 = zeros(ylength-halfrow,xlength-halfcol);
%% Cumulant Martix Calculations
%%%%%%%%%%%%%%%%%%%%%Pull rand line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = halfrow+1:ylength-halfrow;
    for col = halfcol+1:xlength-halfcol;
        
        Var(1:base,1:base) = hhref(row-halfrow:row+halfrow,col-halfcol:col+halfcol);
        sh1(1,1:Window_optimal) = reshape(Var,1,Window_optimal);
        
        Var(1:base,1:base) = hhoff(row-halfrow:row+halfrow,col-halfcol:col+halfcol);
        sh2(1,1:Window_optimal) = reshape(Var,1,Window_optimal);
        
        Var(1:base,1:base) = vvref(row-halfrow:row+halfrow,col-halfcol:col+halfcol);
        sv1(1,1:Window_optimal) = reshape(Var,1,Window_optimal);
        
        Var(1:base,1:base) = vvoff(row-halfrow:row+halfrow,col-halfcol:col+halfcol);
        sv2(1,1:Window_optimal) = reshape(Var,1,Window_optimal);
        
        Var(1:base,1:base) = xxref(row-halfrow:row+halfrow,col-halfcol:col+halfcol);
        sx1(1,1:Window_optimal) = reshape(Var,1,Window_optimal);
        
        Var(1:base,1:base) = xxoff(row-halfrow:row+halfrow,col-halfcol:col+halfcol);
        sx2(1,1:Window_optimal) = reshape(Var,1,Window_optimal);
        
        s1_Noise(1,:) = sh1;
        s1_Noise(2,:) = sv1;
        s1_Noise(3,:) = sx1;
        
        s2_Noise(1,:) = sh2;
        s2_Noise(2,:) = sv2;
        s2_Noise(3,:) = sx2;
        
       index=[1 1
            2 2
            3 3
            1 2
            1 3
            2 3];
        
        S1_2 = [s1_Noise(1,:)
            s1_Noise(2,:)
            s1_Noise(3,:)];
        
        S2_2 = [s2_Noise(1,:)
            s2_Noise(2,:)
            s2_Noise(3,:)];
        
        S1_4 = [s1_Noise(1,:).*s1_Noise(1,:)
            s1_Noise(2,:).*s1_Noise(2,:)
            s1_Noise(3,:).*s1_Noise(3,:)
            s1_Noise(1,:).*s1_Noise(2,:)
            s1_Noise(1,:).*s1_Noise(3,:)
            s1_Noise(2,:).*s1_Noise(3,:)];
        
        S2_4 = [s2_Noise(1,:).*s2_Noise(1,:)
            s2_Noise(2,:).*s2_Noise(2,:)
            s2_Noise(3,:).*s2_Noise(3,:)
            s2_Noise(1,:).*s2_Noise(2,:)
            s2_Noise(1,:).*s2_Noise(3,:)
            s2_Noise(2,:).*s2_Noise(3,:)];
        
        R1_2 = S1_2*S1_2'/Window_optimal;
        R2_2 = S1_2*S2_2'/Window_optimal;
        
        R10_2 = S1_2*S1_2.'/Window_optimal;
        R20_2 = S2_2*S2_2.'/Window_optimal;
        
        R1_4 = S1_4*S1_4'/Window_optimal;
        R2_4 = S1_4*S2_4'/Window_optimal;
        
        for m = 1:6
            for n = 1:6
                
                R1_4(m,n)= R1_4(m,n) - R10_2(index(m,1),index(m,2))*R10_2(index(n,1),index(n,2))...
                    - R1_2(index(m,1),index(n,1))*R1_2(index(m,2),index(n,2))...
                    - R1_2(index(m,1),index(n,2))*R1_2(index(m,2),index(n,1));
                
                R2_4(m,n)= R2_4(m,n) - R10_2(index(m,1),index(m,2))*R20_2(index(n,1),index(n,2))...
                    - R2_2(index(m,1),index(n,1))*R2_2(index(m,2),index(n,2))...
                    - R2_2(index(m,1),index(n,2))*R2_2(index(m,2),index(n,1));
                
            end
        end
        
        [~,eigenvalCov_4] = eig(R1_4,'nobalance');
        
        eye_4 = mean(diag(eigenvalCov_4))/max(abs(diag(eigenvalCov_4)))^2;
        
        [eigenvec_4,eigenval_4] = eig((pinv(R1_4 + eye_4*eye(6)))*R2_4,'nobalance');
        
        [~,srt_4] = sort(abs(diag(eigenval_4)),'descend');
        
        LeigTemp  = (abs(eigenval_4(srt_4(1),srt_4(1))))^2....
            *(abs(eigenvec_4(3,srt_4(1)))^2....
            + abs(eigenvec_4(5,srt_4(1)))^2....
            + abs(eigenvec_4(6,srt_4(1)))^2);
        
        SLeigTemp = (abs(eigenval_4(srt_4(2),srt_4(2))))^2....
            *(abs(eigenvec_4(3,srt_4(2)))^2....
            + abs(eigenvec_4(5,srt_4(2)))^2....
            + abs(eigenvec_4(6,srt_4(2)))^2);
        
        if LeigTemp >= SLeigTemp
            
            vegetation_4(row,col) = eigenval_4(srt_4(1),srt_4(1));
            
            ground_4(row,col) = eigenval_4(srt_4(2),srt_4(2));
            
        else
            
            vegetation_4(row,col) = eigenval_4(srt_4(2),srt_4(2));
            ground_4(row,col) = eigenval_4(srt_4(1),srt_4(1));
        end
    end
end
%% Plotting gv Results
figure(1); imagesc(angle(ground_4)); title('4th Order Ground Interferometric Phase');colorbar;
figure(2); imagesc(angle(vegetation_4)); title('4th Order Vegetation Interferometric Phase');colorbar;
figure(3); imagesc(angle(vegetation_4) - angle(ground_4)); title('4th Order V-G');colorbar;
figure(4); imagesc(abs(ground_4)); title('4th Order Ground Magnitude');colorbar;
figure(5); imagesc(abs(vegetation_4)); title('4th Order Vegetation Magnitude');colorbar;
