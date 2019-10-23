
warning off
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/denoisingsimp_1_noise_models')
warning on

N = 128;
name = 'boat';
M0 = load_image(name,256);
M0 = rescale(crop(M0,N));

n = 1024;
name = 'piece-regular';
f0 = rescale( load_signal(name,n) );

sigma = .1;
M = M0 + randn(N,N)*sigma;
f = f0 + randn(n,1)*sigma;

clf;
subplot(3,1,1);
plot(f0); axis([1 n 0 1]);
title('Clean signal');
subplot(3,1,2);
plot(f-f0); axis([1 n -3*sigma 3*sigma]);
title('Noise');
subplot(3,1,3);
plot(f); axis([1 n 0 1]);
title('Noisy signal');

clf;
imageplot(M0, 'Clean image', 1,3,1);
imageplot(M-M0, 'Noise', 1,3,2);
imageplot(clamp(M), 'Noisy image', 1,3,3);

nbins = 51;
[h,t] = hist( M(:)-M0(:), nbins ); h = h/sum(h);
subplot(3,1,2);
bar(t,h);
axis([-sigma*5 sigma*5 0 max(h)*1.01]);

a = sqrt(3)*sigma;
M = M0 + 2*(rand(N,N)-.5)*a;
f = f0 + 2*(rand(n,1)-.5)*a;

clf;
subplot(3,1,1);
plot(f0); axis([1 n 0 1]);
title('Clean signal');
subplot(3,1,2);
plot(f-f0); axis([1 n -3*sigma 3*sigma]);
title('Noise');
subplot(3,1,3);
plot(f); axis([1 n 0 1]);
title('Noisy signal');

clf;
imageplot(M0, 'Clean image', 1,3,1);
imageplot(M-M0, 'Noise', 1,3,2);
imageplot(clamp(M), 'Noisy image', 1,3,3);

nbins = 51;
[h,t] = hist( M(:)-M0(:), nbins ); h = h/sum(h);
subplot(3,1,2);
bar(t,h);
axis([-sigma*5 sigma*5 0 max(h)*1.01]);

W = log(rand(N,N)).*sign(randn(N,N)); 
W = W/std(W(:))*sigma;
M = M0 + W;

W = log(rand(n,1)).*sign(randn(n,1)); 
W = W/std(W(:))*sigma;
f = f0 + W;

clf;
subplot(3,1,1);
plot(f0); axis([1 n 0 1]);
title('Clean signal');
subplot(3,1,2);
plot(f-f0); axis([1 n -3*sigma 3*sigma]);
title('Noise');
subplot(3,1,3);
plot(f); axis([1 n 0 1]);
title('Noisy signal');

clf;
imageplot(M0, 'Clean image', 1,3,1);
imageplot(M-M0, 'Noise', 1,3,2);
imageplot(clamp(M), 'Noisy image', 1,3,3);

nbins = 51;
[h,t] = hist( M(:)-M0(:), nbins ); h = h/sum(h);
subplot(3,1,2);
bar(t,h);
axis([-sigma*5 sigma*5 0 max(h)*1.01]);

n = 4096;

rho = .05;

x0 = rand(n,1)<rho;

x0 = 2 * x0 .* ( rand(n,1)-.5 );

sigma = .1;
x = x0 + randn(size(x0))*sigma;

clf;
subplot(2,1,1);
plot(x0); axis([1 n -1 1]);
set_graphic_sizes([], 20);
title('Original signal');
subplot(2,1,2);
plot(x); axis([1 n -1 1]);
set_graphic_sizes([], 20);
title('Noisy signal');

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

n = 256;
M0 = rescale(load_image('boat', n));

sigma = 0.06;
M = M0 + randn(n,n)*sigma;

H = M;
H = (H(1:n-1,:) - H(2:n,:))'/sqrt(2);
H = (H(1:n-1,:) - H(2:n,:))'/sqrt(2);

clf;
imageplot(clamp(M), 'Noisy image', 1,2,1);
imageplot(H, 'Derivative image', 1,2,2);

[h,t] = hist(H(:), 100);
h = h/sum(h);

clf;
bar(t, h);
axis([-.5 .5 0 max(h)]);

sigma_est = mad(H(:),1)/0.6745;
disp( strcat(['Estimated noise level=' num2str(sigma_est), ', true=' num2str(sigma)]) );
