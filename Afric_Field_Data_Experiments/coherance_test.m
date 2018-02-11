[hhref,hhoff,vvref,vvoff,xxref,xxoff,P] =  openSARdata();
[xlength, ylength]=size(hhref);
c = zeros([xlength,ylength]);
PolInSAR= zeros([xlength,ylength]);
%%%%%%%%%%%%%%%%%%%%%Pull rand line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = 1:xlength;
    for col = 1:ylength;
        s1 = [hhref(row,col),vvref(row,col),xxref(row,col)];
        s2 = [hhoff(row,col),vvoff(row,col),xxoff(row,col)];
        c(row,col) = (s1*s2')/(sqrt(s1*s1')*sqrt(s2*s2'));
        
        k1 = 1/sqrt(2).*[P.hr(row,col)+P.vr(row,col);P.xr(row,col)-P.vr(row,col);2.*P.xo(row,col)];
        k2 = 1/sqrt(2).*[P.ho(row,col)+P.vo(row,col);P.xo(row,col)-P.vo(row,col);2.*P.xo(row,col)];
        PolInSAR(row,col) = (k1'*k2)/(sqrt(k1'*k1)*sqrt(k2'*k2));
        R1 = k1*k1';
        R2 = k1*k2';
        
%         [U,Uv] = eig(R1*R2');
%         [UU,kk]=sort(abs(angle(diag(Uv))),'ascend');
%         c(row,col) = Uv(kk(1),kk(1));
    end
end
figure(1)
imagesc(angle(PolInSAR));
figure(2)
imagesc(angle(c));