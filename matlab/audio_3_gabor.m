
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/audio_3_gabor')

wlist = 32*[4 8 16 32];  
L = length(wlist);

K = 2;
qlist = wlist/K;

disp( strcat(['Approximate redundancy of the dictionary=' num2str(K*L) '.']) );

n = 1024*32;
options.n = n;
[x0,fs] = load_sound('glockenspiel', n);

options.multichannel = 0;
S = perform_stft(x0,wlist,qlist, options);

exo1()

%% Insert your code here.

plot_spectrogram(S, x0);

x1 = perform_stft(S,wlist,qlist, options);

e = norm(x0-x1)/norm(x0);
disp(strcat(['Reconstruction error (should be 0) = ' num2str(e, 3)]));

sigma = .05;
x = x0 + sigma*randn(size(x0));

S = perform_stft(x,wlist,qlist, options);

T = sigma;
ST = perform_thresholding(S, T, 'soft');

xT = perform_stft(ST,wlist,qlist, options);

err = snr(x0,xT);
clf
plot_spectrogram(ST, xT);
subplot(length(ST)+1,1,1);
title(strcat(['Denoised, SNR=' num2str(err,3), 'dB']));

exo2()

%% Insert your code here.

lambda = .1;
x1 = x;
S1 = perform_stft(x1,wlist,qlist, options);

r = x - x1;
Sr = perform_stft(r, wlist, qlist, options);
S1 = cell_add(S1, Sr);

S1 = perform_thresholding(S1, lambda, 'soft');

x1 = perform_stft(S1,wlist,qlist, options);

exo3()

%% Insert your code here.

e = snr(x0,xbp);
clf
plot_spectrogram(Sbp, xbp);
subplot(length(Sbp)+1,1,1);
title(strcat(['Denoised, SNR=' num2str(e,3), 'dB']));

n = 1024*32;
options.n = n;
s = 3; % number of sound
p = 2; % number of micros
options.subsampling = 1;
x = zeros(n,3);
[x(:,1),fs] = load_sound('bird', n, options);
[x(:,2),fs] = load_sound('male', n, options);
[x(:,3),fs] = load_sound('glockenspiel', n, options);

x = x./repmat(std(x,1), [n 1]);

theta = linspace(0,pi(),s+1); theta(s+1) = [];
theta(1) = .2;
M = [cos(theta); sin(theta)];

y = x*M';

options.multichannel = 1;
S = perform_stft(y, wlist, qlist, options);

y1 = perform_stft(S, wlist, qlist, options);
disp(strcat(['Reconstruction error (should be 0)=' num2str(norm(y-y1,'fro')/norm(y, 'fro')) '.' ]));

lambdaV = .2;

y1 = y;
S1 = S;
niter = 100;
err = [];

for i=1:niter
    % gradient
    r = y - y1;
    Sr = perform_stft(r, wlist, qlist, options);
    S1 = cell_add(S1, Sr);
    % multi-channel thresholding
    %%% BUG HERE %%%%
    % S1 = perform_thresholding(S1, lambdaV, 'soft-multichannel');
    % update the value of lambdaV to match noise
    y1 = perform_stft(S1,wlist,qlist, options);
end

P1 = []; P = [];
for i=1:length(S)
    Si = reshape( S1{i}, [size(S1{i},1)*size(S1{i},2) 2] );
    P1 = cat(1, P1,  Si);
    Si = reshape( S{i}, [size(S{i},1)*size(S{i},2) 2] );
    P = cat(1, P,  Si);
end
P = [real(P);imag(P)];
P1 = [real(P1);imag(P1)];

p = size(P,1);
m = 10000;
sel = randperm(p); sel = sel(1:m);
clf;
subplot(1,2,1);
plot( P(sel,1),P(sel,2), '.' );
title('Tight frame coefficients');
axis([-10 10 -10 10]);
subplot(1,2,2);
plot( P1(sel,1),P1(sel,2), '.' );
title('Basis Pursuit coefficients');
axis([-10 10 -10 10]);

d  = sqrt(sum(P.^2,2));
d1 = sqrt(sum(P1.^2,2));
I = find( d>.2 );
I1 = find( d1>.2 );

Theta  = mod(atan2(P(I,2),P(I,1)), pi());
Theta1 = mod(atan2(P1(I1,2),P1(I1,1)), pi());

nbins = 150;
[h,t] = hist(Theta, nbins);
h = h/sum(h);
[h1,t1] = hist(Theta1, nbins);
h1 = h1/sum(h1);

clf;
subplot(2,1,1);
bar(t,h); axis('tight');
set_graphic_sizes([], 20);
title('Tight frame coefficients');
subplot(2,1,2);
bar(t1,h1); axis('tight');
set_graphic_sizes([], 20);
title('Sparse coefficients');

exo4()

%% Insert your code here.


