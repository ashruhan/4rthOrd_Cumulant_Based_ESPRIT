%% RAYLEIGH Published

variance = 9;
sigma = 1:variance;
samples = 3000;
V = linspace(0,30,samples);
p_v = zeros(variance,samples);
for i = 1:length(sigma)
p_v(i,:) = (V./(2.*pi.*(sigma(i).^2))).*(exp(-(V.^2)./(2*(sigma(i).^2))));
end
figure(1)
for i = 1:length(sigma)
    
    plot(V,p_v(i,:))
    hold on
    
end
hold off
legend('sigma=1','sigma=2','sigma=3','sigma=4','sigma=5','sigma=6','sigma=7','sigma=8','sigma=9')