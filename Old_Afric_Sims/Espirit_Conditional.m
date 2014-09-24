%% Espirit Algorithm conditional Second order
%% Initializations
[hh,vv,xx] =  openSARdata();
%%%%%This is just a work around for right now
hhref = hh.ref;hhoff = hh.off;
vvref = vv.ref;vvoff = vv.off;
xxref = xx.ref;xxoff = xx.off;
%%%%%This is just a work around for right now
r = 3; c = 3;
L = r+c+1;
[ylength,xlength] = size(hhoff);
%% Espirit Martix Calculations
%%%%%%%%%%%%%%%%%%%%%Pull rand line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = 1:ylength;
    for col = 1:xlength;
        if(row<r+1)&&(col<c+1)     % condition 1
            clear('Var','s1','s2','s3','x','y','Lshape')
            
            Var = hhref(1:row+r,1:col+c);
            [x,y] = size(Var); Lreshape = x*y;
            s1(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = hhoff(1:row+r,1:col+c);
            s2(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvref(1:row+r,1:col+c);
            s1(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvoff(1:row+r,1:col+c);
            s2(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxref(1:row+r,1:col+c);
            s1(3,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxoff(1:row+r,1:col+c);
            s2(3,1:Lreshape) = reshape(Var,1,Lreshape);
        elseif(row<r+1)&&(c+1<col)&&(col<xlength-c)    % condition 2
            clear('Var','s1','s2','x','y','Lshape')
            
            Var = hhref(1:row+r,col-c:col+c);
            [x,y] = size(Var); Lreshape = x*y;
            s1(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = hhoff(1:row+r,col-c:col+c);
            s2(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvref(1:row+r,col-c:col+c);
            s1(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvoff(1:row+r,col-c:col+c);
            s2(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxref(1:row+r,col-c:col+c);
            s1(3,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxoff(1:row+r,col-c:col+c);
            s2(3,1:Lreshape) = reshape(Var,1,Lreshape);
        elseif(row<r+1)&&(col>xlength-c)   % condition 3
            clear('Var','s1','s2','s3','x','y','Lshape')
            
            Var = hhref(1:row+r,xlength-((xlength-col)+c):end);
            [x,y] = size(Var); Lreshape = x*y;
            s1(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = hhoff(1:row+r,xlength-((xlength-col)+c):end);
            s2(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvref(1:row+r,xlength-((xlength-col)+c):end);
            s1(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvoff(1:row+r,xlength-((xlength-col)+c):end);
            s2(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxref(1:row+r,xlength-((xlength-col)+c):end);
            s1(3,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxoff(1:row+r,xlength-((xlength-col)+c):end);
            s2(3,1:Lreshape) = reshape(Var,1,Lreshape);            
        elseif(r+1<row)&&(row<ylength-r)&&(col<c+1)   % condition 4
            clear('Var','s1','s2','s3','x','y','Lshape')
            
            Var = hhref(row-r:row+r,1:col+c);
            [x,y] = size(Var); Lreshape = x*y;
            s1(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = hhoff(row-r:row+r,1:col+c);
            s2(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvref(row-r:row+r,1:col+c);
            s1(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvoff(row-r:row+r,1:col+c);
            s2(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxref(row-r:row+r,1:col+c);
            s1(3,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxoff(row-r:row+r,1:col+c);
            s2(3,1:Lreshape) = reshape(Var,1,Lreshape);              
        elseif(r+1<row)&&(row<ylength-r)&&(c+1<col)&&(col<xlength-c)  % condition 5
            clear('Var','s1','s2','s3','x','y','Lshape')
            
            Var = hhref(row-r:row+r,col-c:col+c);
            [x,y] = size(Var); Lreshape = x*y;
            s1(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = hhoff(row-r:row+r,col-c:col+c);
            s2(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvref(row-r:row+r,col-c:col+c);
            s1(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvoff(row-r:row+r,col-c:col+c);
            s2(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxref(row-r:row+r,col-c:col+c);
            s1(3,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxoff(row-r:row+r,col-c:col+c);
            s2(3,1:Lreshape) = reshape(Var,1,Lreshape);
        elseif(r+1<row)&&(row<ylength-r)&&(col>xlength-c)     % condition 6
            clear('Var','s1','s2','s3','x','y','Lshape')
            
            Var = hhref(row-r:row+r,xlength-((xlength-col)+c):end);
            [x,y] = size(Var); Lreshape = x*y;
            s1(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = hhoff(row-r:row+r,xlength-((xlength-col)+c):end);
            s2(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvref(row-r:row+r,xlength-((xlength-col)+c):end);
            s1(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvoff(row-r:row+r,xlength-((xlength-col)+c):end);
            s2(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxref(row-r:row+r,xlength-((xlength-col)+c):end);
            s1(3,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxoff(row-r:row+r,xlength-((xlength-col)+c):end);
            s2(3,1:Lreshape) = reshape(Var,1,Lreshape);            
        elseif(row>ylength-r)&&(col<c+1)    % condition 7
            clear('Var','s1','s2','s3','x','y','Lshape')
            
            Var = hhref(ylength-((ylength-row)+r):end,1:col+c);
            [x,y] = size(Var); Lreshape = x*y;
            s1(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = hhoff(ylength-((ylength-row)+r):end,1:col+c);
            s2(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvref(ylength-((ylength-row)+r):end,1:col+c);
            s1(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvoff(ylength-((ylength-row)+r):end,1:col+c);
            s2(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxref(ylength-((ylength-row)+r):end,1:col+c);
            s1(3,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxoff(ylength-((ylength-row)+r):end,1:col+c);
            s2(3,1:Lreshape) = reshape(Var,1,Lreshape);             
        elseif(row>ylength-r)&&(c+1<col)&&(col<xlength-c)   % condition 8
            clear('Var','s1','s2','s3','x','y','Lshape')
            
            Var = hhref(ylength-((ylength-row)+r):end,col-c:col+c);
            [x,y] = size(Var); Lreshape = x*y;
            s1(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = hhoff(ylength-((ylength-row)+r):end,col-c:col+c);
            s2(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvref(ylength-((ylength-row)+r):end,col-c:col+c);
            s1(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvoff(ylength-((ylength-row)+r):end,col-c:col+c);
            s2(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxref(ylength-((ylength-row)+r):end,col-c:col+c);
            s1(3,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxoff(ylength-((ylength-row)+r):end,col-c:col+c);
            s2(3,1:Lreshape) = reshape(Var,1,Lreshape);                  
        elseif(row>ylength-r)&&(col>xlength-c)   % condition 9
            clear('Var','s1','s2','s3','x','y','Lshape')
            
            Var = hhref(ylength-((ylength-row)+r):end,xlength-((xlength-col)+c):end);
            [x,y] = size(Var); Lreshape = x*y;
            s1(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = hhoff(ylength-((ylength-row)+r):end,xlength-((xlength-col)+c):end);
            s2(1,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvref(ylength-((ylength-row)+r):end,xlength-((xlength-col)+c):end);
            s1(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = vvoff(ylength-((ylength-row)+r):end,xlength-((xlength-col)+c):end);
            s2(2,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxref(ylength-((ylength-row)+r):end,xlength-((xlength-col)+c):end);
            s1(3,1:Lreshape) = reshape(Var,1,Lreshape);
            
            Var = xxoff(ylength-((ylength-row)+r):end,xlength-((xlength-col)+c):end);
            s2(3,1:Lreshape) = reshape(Var,1,Lreshape);              
        end
        
        R1 = s1*s1';  R2 = s1*s2';   A = pinv(R1)*R2;
        
        [~,uv] = eig(A);
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