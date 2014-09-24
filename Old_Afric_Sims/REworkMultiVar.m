[hhref,hhoff,vvref,vvoff,xxref,xxoff] =  openSARdata();

HH = hhref./hhoff;
VV = zeros([10,10]);
s1 = zeros([1,3]);
s2 = zeros([1,3]);
for row = 250:260;
    for col = 400:410;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%To create S1%%%%%%%%%%%%%%%%
        s1(1,1) = vvref(row,col);
        s1(1,2) = xxref(row,col);
        s1(1,3) = hhref(row,col);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%End create S1%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%To create S2%%%%%%%%%%%%%%%%%%%%%%%%%%
        s2(1,1) = vvoff(row,col);
        s2(1,2) = xxoff(row,col);
        s2(1,3) = hhoff(row,col);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End create S1%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        R1 = s1*pinv(s1);
        R2 = s2*pinv(s1);
        [U ,Uv] = eig(R2*pinv(R1));
        %%%%%%%%%%%%%%%%%%%%%%%%%Sorting out Ground from Vege%%%%%%%%%%%%%%%%%%%%
        VV(row-249,col-399) = Uv;
        
    end
end

VVabs = abs(VV)/sqrt(100);
VVangle = angle(VV)/sqrt(100);
[pxabs,pyabs] = gradient(VVabs);
[pxangle,pyangle] = gradient(VVangle);
%%%%%%%%%%%%%%%%%%%%%%%%%%plot%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subplot(2,1,1)
contour(HH)
subplot(2,1,2)
quiver(pxangle,pyangle)


