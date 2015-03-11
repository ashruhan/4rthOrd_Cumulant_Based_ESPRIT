function [hh,vv,xx] =  Copy_of_openSARdata()

fid(1)=fopen('ref_hh.dat','r');fid(2)=fopen('ref_vv.dat','r');
fid(3)=fopen('ref_xx.dat','r');fid(4)=fopen('offset_hh.dat','r');
fid(5)=fopen('offset_vv.dat','r');fid(6)=fopen('offset_xx.dat','r');

h1=fread(fid(1),[730 730],'single');h2=fread(fid(4),[730 730],'single');
v1=fread(fid(2),[730 730],'single');v2=fread(fid(5),[730 730],'single');
x1=fread(fid(3),[730 730],'single');x2=fread(fid(6),[730 730],'single');

for close = 1:6;
fclose(fid(close));
end

hh.ref=h1(1:2:729,:) + 1i*h1(2:2:730,:); hh.off=h2(1:2:729,:) + 1i*h2(2:2:730,:);
vv.ref=v1(1:2:729,:) + 1i*v1(2:2:730,:); vv.off=v2(1:2:729,:) + 1i*v2(2:2:730,:);
xx.ref=x1(1:2:729,:) + 1i*x1(2:2:730,:); xx.off=x2(1:2:729,:) + 1i*x2(2:2:730,:);