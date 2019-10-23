
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/segmentation_3_snakes_levelset')

n = 200; 
[Y,X] = meshgrid(1:n,1:n);

r = n/3;

c = [r r] + 10;

phi1 = sqrt( (X-c(1)).^2 + (Y-c(2)).^2 ) - r;

exo1()

%% Insert your code here.

clf;
subplot(1,2,1);
plot_levelset(phi1);
subplot(1,2,2);
plot_levelset(phi2);

exo2()

%% Insert your code here.

Tmax = 200;

tau = .5;

niter = round(Tmax/tau);

options.order = 2;

phi = phi0;

g0 = grad(phi,options);

d = max(eps, sqrt(sum(g0.^2,3)) );

g = g0 ./ repmat( d, [1 1 2] );

K = -d .* div( g,options );

phi = phi - tau*K;

exo3()

%% Insert your code here.

phi = phi0.^3;

phi1 = perform_redistancing(phi0);

clf;
subplot(1,2,1);
plot_levelset(phi);
title('Before redistancing');
subplot(1,2,2);
plot_levelset(phi1);
title('After redistancing');

n = 200;

name = 'cortex';
f0 = rescale( sum( load_image(name, n), 3) );

g = grad(f0,options);
d0 = sqrt(sum(g.^2,3));

a = 5;

d = perform_blurring( d0,a );

epsilon = 1e-1;

W = 1./(epsilon+d);
W = rescale(-d,.1,1);

clf;
imageplot(f0,'Image to segment',1,2,1);
imageplot(W,'Weight',1,2,2);

exo4()

%% Insert your code here.

clf;
plot_levelset(phi0,0,f0);

tau = .4;

Tmax = 1500;
niter = round(Tmax/tau);

phi = phi0;

gW = grad(W,options);

exo5()

%% Insert your code here.

phi = phi - tau*G;

phi = perform_redistancing(phi);

exo6()

%% Insert your code here.

exo7()

%% Insert your code here.

lambda = 0.8;

c1 = 0.7; 
c2 = 0;

tau = .4;

Tmax = 100;
niter = round(Tmax/tau);

phi = phi0;

exo8()

%% Insert your code here.

phi = phi + tau*G;

exo9()

%% Insert your code here.

exo10()

%% Insert your code here.
