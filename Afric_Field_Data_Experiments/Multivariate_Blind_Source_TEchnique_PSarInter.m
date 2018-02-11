clc;
clear;

[hhref,hhoff,vvref,vvoff,xxref,xxoff] =  openSARdata();

xlength = 365; ylength = 1662;
VV = zeros([xlength,ylength]);
GG = zeros([xlength,ylength]);
NN = zeros([xlength,ylength]);
s1 = zeros([3,49]);
s2 = zeros([3,49]);
 t =1;
%%%%%%%%%%%%%%%%%%%%%Pull rand line%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = fliplr(5:xlength-5);      
      for col = 5:ylength-5; 
          
for d = 1:3;
    h = 1;
    for r = -3:3;
    for c = -3:3;
        if d==1;
        s1(d,h) = hhref(row+r,col+c);
        s2(d,h) = hhoff(row+r,col+c);
         h = h+1;
        
        elseif d==2;
        s1(d,h) = xxref(row+r,col+c);
        s2(d,h) = xxoff(row+r,col+c);
         h = h+1;
        
        elseif d==3;
        s1(d,h) = hhref(row+r,col+c);
        s2(d,h) = hhoff(row+r,col+c);
         h = h+1;
        end
       
    end
    end
        
end
%%%%%%%%%%%%%%%%%%Super Resolution%%%%%%%%%%%%%%%%%%%%%%
 R1 = s1*pinv(s1);
 R2 = s2*pinv(s1);
 [U ,Uv] = eig(R2*pinv(R1));
 %%%%%%%%%%%%%%%%%%%%%%%%%Sorting out Ground from Vege%%%%%%%%%%%%%%%%%%%%
 
 %%%SORTING ROUTING WILL ONLY WORK WHEN THE INTERFEROMIC PHASE DIFFERANCE OF
 %%%GROUND IS LARGER THAN THE PHASE DIFFERANCE OF THE VEGITATION
[UU,kk]=sort(abs(angle(diag(Uv))-pi/2),'ascend');
GG(xlength,ylength) = Uv(kk(1),kk(1));

[UU,kk]=sort(abs(angle(diag(Uv))-pi/4),'ascend');
VV(xlength,ylength) = Uv(kk(1),kk(1));  
      end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%plot%%%%%%%%%%%%%%%%%%%%%%%%%%

