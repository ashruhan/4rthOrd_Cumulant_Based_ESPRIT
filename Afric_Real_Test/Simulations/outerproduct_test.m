N=4;
a=rand(N,1);
tic
 out=bsxfun(@times, a,a');
 out1=bsxfun(@times, a,a.');
toc
%Elapsed time is 0.235199 seconds.
tic;
 out=a*a';
toc
%Elapsed time is 0.325272 seconds.