
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/optim_4b_dr')

set_rand_seeds(123456,123456);

n = 400;

p = round(n/4);

A = randn(p,n) / sqrt(p);

s = 17;

sel = randperm(n);
x0 = zeros(n,1); x0(sel(1:s))=1;

y = A*x0;

proxG = @(x,gamma)max(0,1-gamma./max(1e-15,abs(x))).*x;

t = linspace(-1,1);
plot(t, proxG(t,.3));
axis('equal');

pA = A'*(A*A')^(-1);
proxF = @(x,y)x + pA*(y-A*x);

mu = 1;
gamma = 1;

rproxG = @(x,tau)2*proxG(x,tau)-x;
rproxF = @(x,y)2*proxF(x,y)-x;

niter = 500;

exo1()

%% Insert your code here.

plot(log10(lun(1:end/2)-lun(end)));
axis('tight');

clf;
subplot(2,1,1);
plot_sparse_diracs(x0);
set_graphic_sizes([], 15);
title('Original Signal');
subplot(2,1,2);
plot_sparse_diracs(x);
set_graphic_sizes([], 15);
title('Recovered by L1 minimization');

exo2()

%% Insert your code here.

q = 1000;

slist = 14:2:42;

Slist = slist(mod(0:q-1,length(slist))+1);

U = rand(n,q);
v = sort(U);
v = v( (0:q-1)*n + Slist );
x0 = U<=repmat( v, [n 1] );

y = A*x0;

exo3()

%% Insert your code here.
