function [L] = leverage(N,mu,sigma,T)

clf; % clear the graph figure

t = (0:1:N)'/N; % t is the column vector [0 1/N 2/N ... 1]
W = [0; cumsum(randn(N,1))]/sqrt(N); % S is running sum of N(0,1/N) variables
t = t*T; % adjustment of t
dW = W*sqrt(T); % Brownian Motion dW

S=1*exp(mu*t + sigma*dW); % The GBM stochastic process

Ret=diff(S)./S(1:end-1,:); % The Returns (Ret) of such a process
f=@(x)prod(1+Ret*x); % The product of returns with an argument x
L=fminbnd(@(x)-f(x),-10,10); % Finding the maximum value of x for the product of returns
X=1*(1+L*Ret); % Plugging in the maximum value into the product function
X=cumprod(X); % simulated in blue on the graph figure
S=1*(1+Ret); % The original product of returns for comparison purposes
S=cumprod(S); % simulated in orange on the graph figure

hold on
plot(X);
plot(S);
legend('Optimized Leverage','Original Simulation');
title(['Maximizing Compound Rate under GBM'])
xlabel({['Leverage Factor L = ' num2str(L) ];[' If L>1, then use leverage'];
[' If 1>L>0 then hold the fraction of L'];[' If L<0, then go short. ']})
hold off
end
