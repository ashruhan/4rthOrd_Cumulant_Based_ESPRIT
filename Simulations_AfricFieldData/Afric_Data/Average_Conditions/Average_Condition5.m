%% Condition 5
function [s1,s2] = Average_Condition5(R,C,hh,vv,xx)

Var = hh.ref(R.row-R.r:R.row+R.r,C.col-C.c:C.col+C.c);
[x,y] = size(Var); Lreshape = x*y;
s1.h = reshape(Var,1,Lreshape);

Var = hh.off(R.row-R.r:R.row+R.r,C.col-C.c:C.col+C.c);
s2.h = reshape(Var,1,Lreshape);

Var = vv.ref(R.row-R.r:R.row+R.r,C.col-C.c:C.col+C.c);
s1.v = reshape(Var,1,Lreshape);

Var = vv.off(R.row-R.r:R.row+R.r,C.col-C.c:C.col+C.c);
s2.v = reshape(Var,1,Lreshape);

Var = xx.ref(R.row-R.r:R.row+R.r,C.col-C.c:C.col+C.c);
s1.x = reshape(Var,1,Lreshape);

Var = xx.off(R.row-R.r:R.row+R.r,C.col-C.c:C.col+C.c);
s2.x = reshape(Var,1,Lreshape);

end