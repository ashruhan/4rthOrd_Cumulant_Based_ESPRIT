%% Clearing old values
clc;clear;

%% Init
load('eigenvals_all.mat')

[x,y] = size(one_4);

Lreshape = x*y;

one(1:Lreshape,1) = reshape(one_4,Lreshape,1);
two(1:Lreshape,1) = reshape(two_4,Lreshape,1);
three(1:Lreshape,1) = reshape(three_4,Lreshape,1);
four(1:Lreshape,1) = reshape(four_4,Lreshape,1);
five(1:Lreshape,1) = reshape(five_4,Lreshape,1);
six(1:Lreshape,1) = reshape(six_4,Lreshape,1);
%% Analyzing
figure(1)
subplot(2,3,1)
hist(abs(one))
title('ONE')

% figure(2)
subplot(2,3,2)
hist(abs(two))
title('TWO')

% figure(3)
subplot(2,3,3)
hist(abs(three))
title('THREE')

% figure(4)
subplot(2,3,4)
hist(abs(four))
title('FOUR')

% figure(5)
subplot(2,3,5)
hist(abs(five))
title('FIVE')

% figure(6)
subplot(2,3,6)
hist(abs(six))
title('SIX')