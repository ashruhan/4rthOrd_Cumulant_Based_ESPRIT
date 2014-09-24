%% Condition 9
function [s1,s2] = Average_Condition9(R,C,hh,vv,xx)

Var = hh.ref(R.ylength-((R.ylength-R.row)+R.r):end,C.xlength-((C.xlength-C.col)+C.c):end);
[x,y] = size(Var); Lreshape = x*y;
s1.h = reshape(Var,1,Lreshape);

Var = hh.off(R.ylength-((R.ylength-R.row)+R.r):end,C.xlength-((C.xlength-C.col)+C.c):end);
s2.h = reshape(Var,1,Lreshape);

Var = vv.ref(R.ylength-((R.ylength-R.row)+R.r):end,C.xlength-((C.xlength-C.col)+C.c):end);
s1.v = reshape(Var,1,Lreshape);

Var = vv.off(R.ylength-((R.ylength-R.row)+R.r):end,C.xlength-((C.xlength-C.col)+C.c):end);
s2.v = reshape(Var,1,Lreshape);

Var = xx.ref(R.ylength-((R.ylength-R.row)+R.r):end,C.xlength-((C.xlength-C.col)+C.c):end);
s1.x = reshape(Var,1,Lreshape);

Var = xx.off(R.ylength-((R.ylength-R.row)+R.r):end,C.xlength-((C.xlength-C.col)+C.c):end);
s2.x = reshape(Var,1,Lreshape);
end
