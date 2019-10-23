
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/optim_4_fb')

N = 1024;

s = 5;

t = (-N/2:N/2-1)';
h = (1-t.^2/s^2).*exp( -(t.^2)/(2*s^2) );
h = h-mean(h);

h1 = fftshift(h); % Recenter the filter for fft use.
Phi = @(u)real(ifft(fft(h1).*fft(u)));

hf = real(fftshift(fft(h1))) / sqrt(N);

q = 200;
clf;
subplot(2,1,1);
plot(-N/2+1:N/2, h);
axis([-q q min(h) max(h)]);
title('Filter, Spacial (zoom)');
subplot(2,1,2);
plot(-N/2+1:N/2, hf);
axis([-q q 0 max(hf)]);
title('Filter, Fourier (zoom)');

s = round(N*.03);

rand('state', 1);
randn('state', 1);

sel = randperm(N); sel = sel(1:s);

x0 = zeros(N,1); x0(sel) = 1;
x0 = x0 .* sign(randn(N,1)) .* (1-.3*rand(N,1));

sigma = .06;

y = Phi(x0) + sigma*randn(N,1);

clf; ms = 20;
subplot(2,1,1);
u = x0; u(x0==0) = NaN;
stem(u, 'b.', 'MarkerSize', ms); axis('tight');
title('Signal x_0');
subplot(2,1,2);
plot(y); axis('tight');
title('Measurements y');

lambda = 1;

proxg = @(x,gamma)perform_thresholding(x, lambda*gamma, 'soft');

gradf = @(x)Phi(Phi(x)-y);

L = max(abs(fft(h)))^2;

gamma = 1.95 / L;

x = y;

x = proxg( x - gamma*gradf(x), gamma );

exo1()

%% Insert your code here.

clf;
subplot(2,1,1);
u = x0; u(x0==0) = NaN;
stem(u, 'b.', 'MarkerSize', ms); axis('tight');
title(['Signal x0']);
subplot(2,1,2);
u = x; u(x==0) = NaN;
stem(u, 'b.', 'MarkerSize', ms); axis('tight');

gamma = 1/L;

exo2()

%% Insert your code here.

gamma = 1/L;

exo3()

%% Insert your code here.
