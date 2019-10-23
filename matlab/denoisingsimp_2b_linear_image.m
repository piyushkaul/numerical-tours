
warning off
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/denoisingsimp_2b_linear_image')
warning on

n = 256;
N = n^2;

name = 'hibiscus';
x0 = rescale( sum(load_image(name,n),3) );

clf;
imageplot(x0);

sigma = .08;

y = x0 + sigma*randn(size(x0));

clf;
imageplot(clamp(y));

cconv = @(a,b)real(ifft2( fft2(a).*fft2(b) ));

normalize = @(h)h/sum(h(:));
t = [0:n/2-1, -n/2:-1]';
[Y,X] = meshgrid(t,t);
h = @(mu)normalize( exp( -(X.^2+Y.^2)/(2*mu^2) ) );

mu = 10;
clf;
subplot(2,1,1);
imageplot( fftshift( h(mu) )  ); axis('tight');
title('h');
subplot(2,1,2);
imageplot( fftshift( real(fft2(h(mu))) )  ); axis('tight');
title('fft2(h)');

denoise = @(x,mu)cconv(h(mu), x);

clf;
imageplot( denoise(y,mu) );

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

clf;
imageplot( denoise(y,mu) );

P = 1/N * abs(fft2(x0)).^2;

h_w = real( ifft2( P ./ ( P + sigma^2 ) ) );

clf;
imageplot(crop(fftshift(h_w),n/8)); axis tight;

clf;
imageplot(cconv(y,h_w));
