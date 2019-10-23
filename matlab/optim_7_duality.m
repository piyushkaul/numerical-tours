
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/optim_7_duality')

As = @(u)-div(u);

norm12 = @(u)sum(sum( sqrt(sum( u.^2,3 )) ));
J = @(x)norm12(A(x));

n = 256;

name = 'hibiscus';
x0 = load_image(name,n);
x0 = rescale( sum(x0,3) );

clf;
imageplot(clamp(x0));

sigma = .1;
y = x0 + randn(n,n)*sigma;

clf;
imageplot(clamp(y));

lambda = .2;

mynorm = @(x)norm(x(:));
f = @(x)1/2*mynorm(x-y)^2;
g = @(x)lambda*J(x);
E = @(x)f(x)+g(x);

F = @(u)1/2*mynorm(y-As(u))^2 - 1/2*mynorm(y)^2;

nablaF = @(u)A(As(u)-y);

d = @(u)repmat( sqrt(sum(u.^2,3)), [1 1 2] );
proxG = @(u,gamma)u ./ max( d(u)/lambda, 1 );

gamma = 1/5;

u = zeros(n,n,2);

u = proxG( u - gamma * nablaF(u), gamma );

x = y - As(u);

exo1()

%% Insert your code here.

clf;
imageplot(clamp(x));
