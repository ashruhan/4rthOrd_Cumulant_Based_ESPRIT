%% Clearing old values

[x,y] = size(one_4);

Lreshape = x*y;

one(1:Lreshape,1) = reshape(one_4,Lreshape,1);
two(1:Lreshape,1) = reshape(two_4,Lreshape,1);
three(1:Lreshape,1) = reshape(three_4,Lreshape,1);
four(1:Lreshape,1) = reshape(four_4,Lreshape,1);
five(1:Lreshape,1) = reshape(five_4,Lreshape,1);
six(1:Lreshape,1) = reshape(six_4,Lreshape,1);
%% Analyzing
% plt = linspace(0,60,100);
% figure(1)
% subplot(2,3,1)
% hist(abs(one),plt)
% title('ONE')
% 
% subplot(2,3,2)
% hist(abs(two),plt)
% title('TWO')
% 
% subplot(2,3,3)
% hist(abs(three),plt)
% title('THREE')
% 
% subplot(2,3,4)
% hist(abs(four),plt)
% title('FOUR')
% 
% subplot(2,3,5)
% hist(abs(five),plt)
% title('FIVE')
% 
% subplot(2,3,6)
% hist(abs(six),plt)
% title('SIX')

figure(2)
image(abs(one_4));
title('one');


figure(3)
image(abs(two_4));
title('two');

figure(4)
image(abs(three_4));
title('three');


figure(5)
image(abs(four_4));
title('four');

figure(6)
image(abs(five_4));
title('five');

figure(7)
image(abs(six_4));
title('six');

figure(8)
image(angle(one_4));
title('one');


figure(9)
image(angle(two_4));
title('two');

figure(10)
image(angle(three_4));
title('three');


figure(11)
image(angle(four_4));
title('four');

figure(12)
image(angle(five_4));
title('five');

figure(13)
image(angle(six_4));
title('six');