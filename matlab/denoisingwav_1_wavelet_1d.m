
warning off
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/denoisingwav_1_wavelet_1d')
warning on

N = 2048*2;

name = 'piece-regular';
f0 = load_signal(name, N);
f0 = rescale(f0,.05,.95);

sigma = 0.05;

f = f0 + randn(size(f0))*sigma;

clf;
subplot(2,1,1);
plot(f0); axis([1 N 0 1]);
title('Clean signal');
subplot(2,1,2);
plot(f); axis([1 N 0 1]);
title('Noisy signal');

Theta0 = @(x,T)x .* (abs(x)>T);

Theta1 = @(x,T)max(0, 1-T./max(abs(x),1e-9)) .* x;

t = linspace(-3,3,1024)'; T = 1;
clf;
plot( t, [Theta0(t,T), Theta1(t,T)], 'LineWidth', 2 );
axis('equal'); axis('tight');
legend('\Theta^0', '\Theta^1');

options.ti = 0; Jmin = 4;
W  = @(f) perform_wavelet_transf(f,Jmin,+1,options);
Wi = @(fw)perform_wavelet_transf(fw,Jmin,-1,options);

x = W(f);

x1 = Theta0(x, 3*sigma);

clf;
subplot(2,1,1);
plot_wavelet(x,Jmin); axis([1 N -1 1]);
title('W(f)');
subplot(2,1,2);
plot_wavelet(Theta0(W(f),T),Jmin); axis([1 N -1 1]);
title('\Theta_0(W(f))');

f1 = Wi(x1);

clf;
subplot(2,1,1);
plot(f); axis([1 N 0 1]);
title('f');
subplot(2,1,2);
plot(f1); axis([1 N 0 1]);
title('f_1');

x = W(f); 
reinject = @(x1)assign(x1, 1:2^Jmin, x(1:2^Jmin));

Theta0W = @(f,T)Wi(Theta0(W(f),T));
Theta1W = @(f,T)Wi(reinject(Theta1(W(f),T)));

exo1()

%% Insert your code here.

clf;
subplot(2,1,1);
plot(fBest0); axis([1 N 0 1]); e = snr(fBest0,f0);
title(['Hard, SNR=' num2str(e,3) 'dB']);
subplot(2,1,2);
plot(fBest1); axis([1 N 0 1]); e = snr(fBest1,f0);
title(['Soft, SNR=' num2str(e,3) 'dB']);

options.ti = 1;
W  = @(f) perform_wavelet_transf(f,Jmin,+1,options);
Wi = @(fw)perform_wavelet_transf(fw,Jmin,-1,options);

fw = W(f);

nJ = size(fw,3)-4;
clf;
subplot(5,1, 1);
plot(f0); axis('tight');
title('Signal');
i = 0;
for j=1:3
    i = i+1;
    subplot(5,1,i+1);
    plot(fw(:,1,nJ-i+1)); axis('tight');
    title(strcat(['Scale=' num2str(j)]));
end
subplot(5,1, 5);
plot(fw(:,1,1)); axis('tight');
title('Low scale');

x = W(f); 
reinject = @(x1)assign(x1, 1:N, x(1:N));

Theta0W = @(f,T)Wi(Theta0(W(f),T));
Theta1W = @(f,T)Wi(reinject(Theta1(W(f),T)));

exo2()

%% Insert your code here.

clf;
subplot(2,1,1);
plot(fBest0); axis([1 N 0 1]); e = snr(fBest0,f0);
title(['Hard TI, SNR=' num2str(e,3) 'dB']);
subplot(2,1,2);
plot(fBest1); axis([1 N 0 1]); e = snr(fBest1,f0);
title(['Soft TI, SNR=' num2str(e,3) 'dB']);
