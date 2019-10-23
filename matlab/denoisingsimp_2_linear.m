
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/denoisingsimp_2_linear')

N = 1024;

name = 'piece-regular';
x0 = rescale( load_signal(name,N) );

sigma = .04;

y = x0 + sigma*randn(size(x0));

clf;
subplot(2,1,1);
plot(x0); axis([1 N -.05 1.05]);
subplot(2,1,2);
plot(y); axis([1 N -.05 1.05]);

cconv = @(a,b)real(ifft( fft(a).*fft(b) ));

normalize = @(h)h/sum(h(:));
t = [0:N/2-1, -N/2:-1]';
h = @(mu)normalize( exp( -(t.^2)/(2*mu^2) ) );

mu = 10;
clf;
subplot(2,1,1);
plot( t, h(mu)  ); axis('tight');
title('h');
subplot(2,1,2);
plot( t, real(fft(h(mu)))  ); axis('tight');
title('fft(h)');

denoise = @(x,mu)cconv(h(mu), x);

clf;
plot( denoise(y,mu) );
axis([1 N -.05 1.05]);

exo1()

%% Insert code here

exo2()

%% Insert code here

clf;
plot( denoise(y,mu) );
axis([1 N -.05 1.05]);

P = 1/N * abs(fft(x0)).^2;

h_w = real( ifft( P ./ ( P + sigma^2 ) ) );

clf;
plot(fftshift(h_w)); axis tight;

clf;
plot(cconv(y,h_w));
axis([1 N -.05 1.05]);
