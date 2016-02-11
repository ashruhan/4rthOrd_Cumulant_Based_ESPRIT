%% Espirit Algorithm conditional Second order
%% Initializations
[hh,vv,xx] =  openSARdata();
r = 3; c = 3;

[ylength,xlength] = size(hh.off);
g = zeros(ylength,xlength);
v = zeros(ylength,xlength);
n = zeros(ylength,xlength);

%% Espirit Martix Calculations
%%%%%%%%%%%%%%%%%%%%%Pull rand line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = 1:ylength;
    for col = 1:xlength;
        
        clear('R1','R2','S1','S2','A','uv','kk')
        R.r = r;R.row = row;R.ylength = ylength;
        C.c = c;C.col = col;C.xlength = xlength;
        
        if(row<r+1)&&(col<c+1)     % condition 1
            [s1,s2] = Average_Condition1(R,C,hh,vv,xx);
            
        elseif(row<r+1)&&(c+1<col)&&(col<xlength-c)    % condition 2
            [s1,s2] = Average_Condition2(R,C,hh,vv,xx);
            
        elseif(row<r+1)&&(col>xlength-c)   % condition 3
            [s1,s2] = Average_Condition3(R,C,hh,vv,xx);
            
        elseif(r+1<row)&&(row<ylength-r)&&(col<c+1)   % condition 4
            [s1,s2] = Average_Condition4(R,C,hh,vv,xx);
            
        elseif(r+1<row)&&(row<ylength-r)&&(c+1<col)&&(col<xlength-c)  % condition 5
            [s1,s2] = Average_Condition5(R,C,hh,vv,xx);
            
        elseif(r+1<row)&&(row<ylength-r)&&(col>xlength-c)     % condition 6
            [s1,s2] = Average_Condition6(R,C,hh,vv,xx);
            
        elseif(row>ylength-r)&&(col<c+1)    % condition 7
            [s1,s2] = Average_Condition7(R,C,hh,vv,xx);
            
        elseif(row>ylength-r)&&(c+1<col)&&(col<xlength-c)   % condition 8
            [s1,s2] = Average_Condition8(R,C,hh,vv,xx);
            
        elseif(row>ylength-r)&&(col>xlength-c)   % condition 9
            [s1,s2] = Average_Condition9(R,C,hh,vv,xx);
        end
        
        S1(1,:) = s1.h; S2(1,:) = s2.h;
        S1(2,:) = s1.v; S2(2,:) = s2.v;
        S1(3,:) = s1.x; S2(3,:) = s2.x;
        
        
        R1 = S1*S1';  R2 = S1*S2';
        
        [eigenvector,eigenvalue] = eig(pinv(R1)*R2);
        
        [~,kk]=sort(angle(diag(eigenvalue)),'descend');
        
        
        if (abs(eigenvector(1,kk(1)))^2+abs(eigenvector(2,kk(1)))^2)<...
                (abs(eigenvector(1,kk(2)))^2+abs(eigenvector(2,kk(2)))^2)
            
            g(row,col) = eigenvalue(kk(1),kk(1));
            v(row,col) = eigenvalue(kk(2),kk(2));
        else
            
            g(row,col) = eigenvalue(kk(2),kk(2));
            v(row,col) = eigenvalue(kk(1),kk(1));
            
        end
               
        n(row,col) = eigenvalue(kk(3),kk(3));
    end
end

%% Plotting Results
figure(4); imagesc(angle(g)); title('angle(g)');
figure(5); imagesc(angle(v)); title('angle(v)');
figure(6); imagesc(angle(n)); title('angle(n)');
figure(7); imagesc(abs(g)); title('abs(g)');
figure(8); imagesc(abs(v)); title('abs(v)');
figure(9); imagesc(abs(n)); title('abs(n)');