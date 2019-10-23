
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/segmentation_4_convex_segmentation')

name = 'hibiscus';
n = 256;
I = rescale( load_image(name,n) );

clf;
imageplot(I);

c0 = [1;0;0];
c1 = [0;1;0];

exo1()

%% Insert your code here.

w = w0-w1;

options.bound = 'per';
Grad = @(x)grad(x,options);
Div = @(x)div(x,options);

[X Y] = meshgrid(0:n-1, 0:n-1);
K = 1 + 4*sin(X*pi/n).^2 + 4*sin(Y*pi/n).^2;

Replicate = @(z)deal(z, Grad(z));
ProjC = @(f,u)Replicate( real( ifft2( fft2( f - Div(u) ) ./ K ) ) );

ProxG = @(f,u,gamma)ProjC(f,u);

lambda = .1;

ProxF0 = @(f,gamma)max(0, min(1, f-gamma*w)  );

amplitude = @(u)repmat( sqrt( sum(u.^2, 3) ), [1 1 2]);
ProxL1 = @(u,gamma)max(0,1-gamma./max(1e-9, amplitude(u))) .* u;

ProxF = @(f,u,gamma)deal( ProxF0(f,gamma), ProxL1(u,gamma*lambda) );

mu = 1;
gamma = 1;

niter = 800;

exo2()

%% Insert your code here.

clf;
imageplot(f);

exo3()

%% Insert your code here.
