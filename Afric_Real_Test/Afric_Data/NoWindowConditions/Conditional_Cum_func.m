%% Cumulant Algorithm conditional Fourth order
%% Initializations
[hh,vv,xx] =  openSARdata();
r = 3; c = 3;
[ylength,xlength] = size(hh.off);
%% Cumulant Martix Calculations
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
        
        S1 = [s1.h.*s1.v;
            s1.h.*s1.x;
            s1.v.*s1.x];
        
        S2 = [s2.h.*s2.v;
            s2.h.*s2.x;
            s2.v.*s2.x];
        
        
        R1 = S1*S1';  
        R2 = S1*S2';   
               
        [~,uv] = eig(pinv(R1)*R2);
        
        [~,kk]=sort(angle(diag(uv)),'ascend');
        
        g(row,col) = uv(kk(1),kk(1)); %#ok<*SAGROW>
        v(row,col) = uv(kk(2),kk(2));
        n(row,col) = uv(kk(3),kk(3));
        
    end
end

%% Plotting Results
[x,y] = size(g);

figure(1); histg = reshape(g,x*y,1); hist(angle(histg),100);
figure(2); histv = reshape(v,x*y,1); hist(angle(histv),100)
figure(3); histn = reshape(n,x*y,1); hist(angle(histn),100);

figure(4); imagesc(angle(g)); title('angle(g)');
figure(5); imagesc(angle(v)); title('angle(v)');
figure(6); imagesc(angle(n)); title('angle(n)');
figure(7); imagesc(abs(g)); title('abs(g)');
figure(8); imagesc(abs(v)); title('abs(v)');
figure(9); imagesc(abs(n)); title('abs(n)');