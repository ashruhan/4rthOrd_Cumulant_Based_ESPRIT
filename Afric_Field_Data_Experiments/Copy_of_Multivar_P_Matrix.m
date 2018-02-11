[hhref,hhoff,vvref,vvoff,xxref,xxoff] =  openSARdata();
r = 3; c = 3;
L = r+c+1;
Lreshape = (r+c+1)^2;
xlength = 365; ylength = 1662;
%%%%%%%%%%%%%%%%%%%%%Pull rand line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = r+1:xlength-r;
    for col = c+1:ylength-c;
        
        Var(1:L,1:L) = hhref(row-r:row+r,col-c:col+c); s1(1,1:Lreshape) = reshape(Var,1,Lreshape);
        Var(1:L,1:L) = hhoff(row-r:row+r,col-c:col+c); s2(1,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var(1:L,1:L) = vvref(row-r:row+r,col-c:col+c); s1(2,1:Lreshape) = reshape(Var,1,Lreshape);
        Var(1:L,1:L) = vvoff(row+r,col+c); s2(2,1:Lreshape) = reshape(Var,1,Lreshape);
        
        Var(1:L,1:L) = xxref(row-r:row+r,col-c:col+c);s1(3,1:Lreshape) = reshape(Var,1,Lreshape);
        Var(1:L,1:L) = xxoff(row-r:row+r,col-c:col+c);s2(3,1:Lreshape) = reshape(Var,1,Lreshape);
        
        
        R1 = s1*s1';
        R2 = s1*s2';
        A=R2*R1^(-1);
        [~,uv] = eig(A);
        
        [~,kk]=sort(angle(diag(uv)),'ascend');       

        G.g(row,col) = uv(kk(1),kk(1));
        V.v(row,col) = uv(kk(2),kk(2));
        N.n(row,col) = uv(kk(3),kk(3));   

    end
end



