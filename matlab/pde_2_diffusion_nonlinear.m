
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/pde_2_diffusion_nonlinear')

n = 256;

name = 'hibiscus';
f0 = load_image(name,n);
f0 = rescale( sum(f0,3) );

clf;
imageplot(f0);

cconv = @(f,h)real(ifft2(fft2(f).*fft2(h)));

t = [0:n/2 -n/2+1:-1];
[X2,X1] = meshgrid(t,t);
normalize = @(h)h/sum(h(:));
h = @(sigma)normalize( exp( -(X1.^2+X2.^2)/(2*sigma^2) ) );

sigma = .5;

blur = @(f)cconv(f,h(sigma));

g = @(s,lambda)1./sqrt( 1+(s/lambda).^2 );

amplitude = @(u)repmat( sqrt( sum(u.^2,3) ), [1 1 2]);
A = @(f)amplitude(grad(blur(f)));

f = f0;

lambda = .01;

tau = .2;

f = f + tau * div( g(A(f),lambda) .* grad(f) );

T = .5/lambda;

niter = ceil(T/tau);

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

epsilon = 1e-6;
amplitude = @(u)sqrt(sum(u.^2,3)+epsilon^2);
normalize = @(u)u./repmat( amplitude(u), [1 1 2]);
curv = @(f)div( normalize(grad(f)) );

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.
