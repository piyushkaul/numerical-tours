
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/pde_1_diffusion_linear')

n = 256;

name = 'hibiscus';
f0 = load_image(name,n);
f0 = rescale( sum(f0,3) );

clf;
imageplot(f0);

h = 1/n;
delta = @(f)1/h^2 * div(grad(f));

tau = .5 * h^2/4;

T = 1e-3;

niter = ceil(T/tau);

f = f0;

f = f + tau * delta(f);

exo1()

%% Insert your code here.

cconv = @(f,h)real(ifft2(fft2(f).*fft2(h)));

t = [0:n/2 -n/2+1:-1];
[X2,X1] = meshgrid(t,t);
normalize = @(h)h/sum(h(:));
h = @(t)normalize( exp( -(X1.^2+X2.^2)/(4*t) ) );

heat = @(f, t)cconv(f,h(t));

clf;
imageplot(heat(f0,2));

exo2()

%% Insert your code here.
