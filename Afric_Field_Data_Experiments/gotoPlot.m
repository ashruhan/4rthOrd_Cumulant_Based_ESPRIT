function gotoPlot()
g = evalin('base','g');
v = evalin('base','v');
n = evalin('base','n');

RI = imref2d(size(angle(n)));
RI.XWorldLimits = [2 16];
RI.YWorldLimits = [2 8];

figure(1)
imshow(angle(g),RI);title('angle(g)');
figure(2)
imshow(angle(v),RI);title('angle(v)');
figure(3)
imshow(angle(n),RI);title('angle(n)');
figure(4)
imshow(abs(g),RI);title('abs(g)');
figure(5)
imshow(abs(v),RI);title('abs(v)');
figure(6)
imshow(abs(n),RI);title('abs(n)');