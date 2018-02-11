%% Algorithm Fourth order
%% Initializations

[hh,vv,xx] =  openSARdata();
%%%%%This is just a work around for right now
hhref = hh.ref;hhoff = hh.off;
vvref = vv.ref;vvoff = vv.off;
xxref = xx.ref;xxoff = xx.off;
%%%%%This is just a work around for right now

[ylength,xlength] = size(hhoff);
r = 50; c = 50;
L = r+c+1;
Lreshape = (r+c+1)^2;

% g = zeros(xlength-r,ylength-c);
% v = zeros(xlength-r,ylength-c);
% n = zeros(xlength-r,ylength-c);

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
        
        S1 = [sh1.*sv1;sh1.*sx1;sv1.*sx1];
        S2 = [sh2.*sv2;sh2.*sx2;sv2.*sx2];
        
        
        R1 = S1*S1';  R2 = S1*S2';   A = pinv(R1)*R2;
   
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