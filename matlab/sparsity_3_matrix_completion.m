
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/sparsity_3_matrix_completion')

n = 100;

r = 10;

x0 = randn(n,r)*randn(r,n);

plot(svd(x0), '.-');
axis tight;

P = round( n*log(n)*r*1 );

I = randperm(n*n); I = I(1:P); I = I(:);

Phi  = @(x)x(I);
PhiS= @(y)reshape( accumarray(I, y, [n*n 1], @sum), [n n]);

y = Phi(x0);

ProxF = @(x,gamma)x + PhiS(y-Phi(x));

SoftThresh = @(x,gamma)max(0,1-gamma./max(abs(x),1e-10)).*x;

t = linspace(-10,10,1000);
h = plot(t, SoftThresh(t,3)); axis tight; axis equal;
set(h, 'LineWidth', 2);

prod = @(a,b,c)a*b*c;
SoftThreshDiag = @(a,b,c,gamma)a*diag(SoftThresh(diag(b),gamma))*c';
ProxG = @(x,gamma)apply_multiple_ouput(@(a,b,c)SoftThreshDiag(a,b,c,gamma), @svd, x);

rProxF = @(x,gamma)2*ProxF(x,gamma)-x;
rProxG = @(x,gamma)2*ProxG(x,gamma)-x;

mu = 1;
gamma = 1;

exo1()

%% Insert your code here.

disp(['|A-A_0|/|A_0| = ' num2str(norm(x-x0)/norm(x), 2)]);

exo2()

%% Insert your code here.

alpha = 1;
[U,R] = qr(randn(n));
[V,R] = qr(randn(n));
S = (1:n).^(-alpha);
x0 = U*diag(S)*V';

clf;
h = plot(S); axis tight;
set(h, 'LineWidth', 2);

P = n*n/4;

I = randperm(n*n); I = I(1:P); I = I(:);
Phi  = @(x)x(I);
PhiS= @(y)reshape( accumarray(I, y, [n*n 1], @sum), [n n]);

sigma = std(x0(:))/5;

y = Phi(x0)+sigma*randn(P,1);

lambda = .01;

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.
