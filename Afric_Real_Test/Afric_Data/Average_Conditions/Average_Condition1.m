%% Condition 1
function [s1,s2] = Average_Condition1(R,C,hh,vv,xx)

Var = hh.ref(1:R.row+R.r,1:C.col+C.c);
[x,y] = size(Var); Lreshape = x*y;
s1.h = reshape(Var,1,Lreshape);

Var = hh.off(1:R.row+R.r,1:C.col+C.c);
s2.h = reshape(Var,1,Lreshape);

Var = vv.ref(1:R.row+R.r,1:C.col+C.c);
s1.v = reshape(Var,1,Lreshape);

Var = vv.off(1:R.row+R.r,1:C.col+C.c);
s2.v = reshape(Var,1,Lreshape);

Var = xx.ref(1:R.row+R.r,1:C.col+C.c);
s1.x = reshape(Var,1,Lreshape);

Var = xx.off(1:R.row+R.r,1:C.col+C.c);
s2.x = reshape(Var,1,Lreshape);

end