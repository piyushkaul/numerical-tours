
n = 1024*16;
options.n = n;
[x,fs] = load_sound('bird', n);

sound(x(:)',fs);

clf;
plot(1:n,x);
axis('tight');
set_graphic_sizes([], 20);
title('Signal');

p = 512;
t = 1:n;
clf;
sel = n/4 + (0:p-1);
subplot(2,1,1);
plot(t(sel),x(sel)); axis tight;
sel = n/2 + (0:p-1);
subplot(2,1,2);
plot(t(sel),x(sel)); axis tight;

exo1()

%% Insert your code here.

t = linspace(-10,10,2048);
eta = 1e-5;
vmin = -2;

h = double( abs(t)<1 );
hf = fftshift(abs(fft(h)));
hf = log10(eta+hf); hf = hf/max(hf);
clf;
subplot(2,1,1);
title('Block window');
plot(t, h); axis([-2 2, -.1, 1.1]);
subplot(2,1,2);
plot(t, hf); axis([-2 2, vmin, 1.1]);
title('Fourier transform');

h = cos(t*pi()/2) .* double(abs(t)<1);
hf = fftshift(abs(fft(h)));
hf = log10(eta+hf); hf = hf/max(hf);
clf;
subplot(2,1,1);
title('Hamming window');
plot(t, h); axis([-2 2, -.1, 1.1]);
subplot(2,1,2);
plot(t, hf); axis([-2 2, vmin, 1.1]);
title('Fourier transform');

h = (cos(t*pi())+1)/2 .* double(abs(t)<1);
hf = fftshift(abs(fft(h)));
hf = log10(eta+hf); hf = hf/max(hf);
clf;
subplot(2,1,1);
title('Haning window');
plot(t, h); axis([-2 2, -.1, 1.1]);
subplot(2,1,2);
plot(t, hf); axis([-2 2, vmin, 1.1]);
title('Fourier transform');

h = sqrt(2)/2 * (1+cos(t*pi())) ./ sqrt( 1+cos(t*pi()).^2 ) .* double(abs(t)<1);
hf = fftshift(abs(fft(h)));
hf = log10(eta+hf); hf = hf/max(hf);
clf;
subplot(2,1,1);
title('Normalized Haning window');
plot(t, h); axis([-2 2, -.1, 1.1]);
subplot(2,1,2);
plot(t, hf); axis([-2 2, vmin, 1.1]);
title('Fourier transform');

w = 64*2;

q = w/2;

t = 0:3*w-1;
t1 = t-2*w;
f = w/8;

g1 = sin( pi*t/w ).^2 .* double(t<w);

g2 = sin( pi*t1/w ).^2 .* double( t1<w & t1>=0 );

g3 = g1 .* sin( t * 2*pi/w * f);

g4 = g2 .* sin( t * 2*pi/w * f);

clf;
subplot(2,2,1);
plot(g1); axis('tight');
title('Position 0, frequency 0');
subplot(2,2,2);
plot(g2); axis('tight');
title('Position 2*w, frequency 0');
subplot(2,2,3);
plot(g3); axis('tight');
title('Position 0, frequency w/8');
subplot(2,2,4);
plot(g4); axis('tight');
title('Position 2*w, frequency w/8');

S = perform_stft(x,w,q, options);

clf; imageplot(abs(S)); axis('on');

plot_spectrogram(S,x);

e = norm(x,'fro').^2;

eS = norm(abs(S),'fro').^2;
disp(strcat(['Energy conservation (should be 1)=' num2str(e/eS)]));

x1 = perform_stft(S,w,q, options);
disp(strcat(['Reconstruction error (should be 0)=' num2str( norm(x-x1, 'fro')./norm(x,'fro') ) ]));

sigma = .2;
xn = x + randn(size(x))*sigma;

sound(xn,fs);

clf;
subplot(2,1,1);
plot(x); axis([1 n -1.2 1.2]);
set_graphic_sizes([], 20);
title('Original signal');
subplot(2,1,2);
plot(xn); axis([1 n -1.2 1.2]);
set_graphic_sizes([], 20);
title('Noisy signal');

Sn = perform_stft(xn,w,q, options);
SnT = perform_thresholding(Sn, 2*sigma, 'hard');

subplot(2,1,1);
plot_spectrogram(Sn);
subplot(2,1,2);
plot_spectrogram(SnT);

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

Sn = perform_stft(xn,w,q, options);
SnT = perform_thresholding(Sn, sigma, 'block');

subplot(2,1,1);
plot_spectrogram(Sn);
subplot(2,1,2);
plot_spectrogram(SnT);

exo4()

%% Insert your code here.
