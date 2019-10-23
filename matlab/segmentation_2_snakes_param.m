
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/segmentation_2_snakes_param')

gamma0 = [0.78 0.14 0.42 0.18 0.32 0.16 0.75 0.83 0.57 0.68 0.46 0.40 0.72 0.79 0.91 0.90]' + ...
     1i* [0.87 0.82 0.75 0.63 0.34 0.17 0.08 0.46 0.50 0.25 0.27 0.57 0.73 0.57 0.75 0.79]';

p = 256;

curvabs = @(gamma)[0;cumsum( 1e-5 + abs(gamma(1:end-1)-gamma(2:end)) )];
resample1 = @(gamma,d)interp1(d/d(end),gamma,(0:p-1)'/p, 'linear');
resample = @(gamma)resample1( [gamma;gamma(1)], curvabs( [gamma;gamma(1)] ) );

gamma1 = resample(gamma0);

clf;
h = plot(gamma1([1:end 1]), 'k');
set(h, 'LineWidth', 2); axis('tight'); axis('off');

BwdDiff = @(c)c - c([end 1:end-1]);
FwdDiff = @(c)c([2:end 1]) - c;
dotp = @(c1,c2)real(c1.*conj(c2));

normalize = @(v)v./max(abs(v),eps);
tangent = @(gamma)normalize( FwdDiff(gamma) );
normal = @(gamma)-1i*tangent(gamma);

delta = .03;
gamma2 = gamma1 + delta * normal(gamma1);
gamma3 = gamma1 - delta * normal(gamma1);

clf;
hold on;
h = plot(gamma1([1:end 1]), 'k'); set(h, 'LineWidth', 2); 
h = plot(gamma2([1:end 1]), 'r--'); set(h, 'LineWidth', 2); 
h = plot(gamma3([1:end 1]), 'b--'); set(h, 'LineWidth', 2); 
axis('tight'); axis('off');

normalC = @(gamma)BwdDiff(tangent(gamma)) ./ abs( FwdDiff(gamma) );

dt = 0.001 / 100;

Tmax = 3 / 100;
niter = round(Tmax/dt);

gamma = gamma1;

gamma = gamma + dt * normalC(gamma);

gamma = resample(gamma);

exo1()

%% Insert your code here.

n = 200;
nbumps = 40;
theta = rand(nbumps,1)*2*pi;
r = .6*n/2; a = [.62*n .6*n];
x = round( a(1) + r*cos(theta) );
y = round( a(2) + r*sin(theta) );
W = zeros(n); W( x + (y-1)*n ) = 1;
W = perform_blurring(W,10);
W = rescale( -min(W,.05), .3,1);

clf;
imageplot(W);

options.order = 2;
G = grad(W, options);
G = G(:,:,1) + 1i*G(:,:,2);

EvalG = @(gamma)interp2(1:n,1:n, G, imag(gamma), real(gamma));
EvalW = @(gamma)interp2(1:n,1:n, W, imag(gamma), real(gamma));

r = .98*n/2;
p = 128; % number of points on the curve
theta = linspace(0,2*pi,p+1)'; theta(end) = [];
gamma0 = n/2*(1+1i) +  r*(cos(theta) + 1i*sin(theta));

gamma = gamma0;

dt = 1;

Tmax = 5000;
niter = round(Tmax/dt);

lw = 2;
clf; hold on;
imageplot(W);
h = plot(imag(gamma([1:end 1])),real(gamma([1:end 1])), 'r');
set(h, 'LineWidth', lw);
axis('ij');

N = normal(gamma);
g = - EvalW(gamma).*normalC(gamma) + dotp(EvalG(gamma), N) .* N;
gamma = gamma - dt*g;

gamma = resample( gamma );

exo2()

%% Insert your code here.

n = 256;
f = rescale( sum(load_image('cortex', n), 3 ) );

clf;
imageplot(f);

options.order = 2;
G = grad(f,options);
d = sqrt(sum(G.^2,3));

a = 3;
d = perform_blurring(d,a);

d = min(d,.4);
W = rescale(-d,.8,1);

clf;
imageplot(W);

p = 128;

exo3()

%% Insert your code here.

dt = 2;

Tmax = 9000;
niter = round(Tmax/dt);

exo4()

%% Insert your code here.

n = 256;
f = rescale( sum(load_image('cortex', n), 3 ) );
f = f(46:105,61:120);
n = size(f,1);

clf;
imageplot(f);

exo5()

%% Insert your code here.

x0 = 4 + 55i;
x1 = 53 + 4i;

p = 128;
t = linspace(0,1,p)';
gamma0 = t*x1 + (1-t)*x0;

gamma = gamma0;

clf; hold on;
imageplot(W);
h = plot(imag(gamma([1:end])),real(gamma([1:end])), 'r'); set(h, 'LineWidth', 2);
h = plot(imag(gamma([1 end])),real(gamma([1 end])), 'b.'); set(h, 'MarkerSize', 30);
axis('ij');

curvabs = @(gamma)[0;cumsum( 1e-5 + abs(gamma(1:end-1)-gamma(2:end)) )];
resample1 = @(gamma,d)interp1(d/d(end),gamma,(0:p-1)'/(p-1), 'linear');
resample = @(gamma)resample1( gamma, curvabs(gamma) );

dt = 1/10;

Tmax = 2000*4/7;
niter = round(Tmax/dt);

exo6()

%% Insert your code here.
