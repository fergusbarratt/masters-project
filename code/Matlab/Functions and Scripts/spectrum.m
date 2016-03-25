W = 80:0.2:120;
gamma = 1;
E = 5;
g = 10;
kappa = 10;
W0 = 100;

Y1 = (abs(E) ./ g)^4 * (((2 * (0.5 * (kappa + gamma ./ 2))^3 ./ pi) ./ ((0.5 * (kappa+(gamma ./ 2)))^2 + (W - W0 + g) .^ 2).^2) + 2 * (0.5 * ((kappa+(gamma./2)))^3./pi)./((0.5 * (kappa + gamma ./ 2))^2 + (W-W0-g).^2).^2);

Y2 = 0.5*((0.5*(kappa+gamma./2)./pi)./(0.25*(kappa+gamma/2)^2+(W-W0+g).^2) + (0.5*(kappa+gamma/2)/pi)./(0.25*(kappa+gamma/2)^2+(W - W0 - g).^2));

plot(W, Y1)
hold on
plot(W, Y2)
hold off
