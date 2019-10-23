
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/audio_2_separation')

extend_stack_size(4);

n = 1024*16;
s = 3; % number of sound
p = 2; % number of micros
options.subsampling = 1;
x = zeros(n,3);
[x(:,1),fs] = load_sound('bird', n, options);
[x(:,2),fs] = load_sound('female', n, options);
[x(:,3),fs] = load_sound('male', n, options);

x = x./repmat(std(x,1), [n 1]);

theta = linspace(0,pi(),s+1); theta(s+1) = [];
theta(1) = .2;
M = [cos(theta); sin(theta)];

y = x*M';

clf;
for i=1:s
    subplot(s,1,i);
    plot(x(:,i)); axis('tight');
    set_graphic_sizes([], 20);
    title(strcat('Source #',num2str(i)));
end

clf;
for i=1:p
    subplot(p,1,i);
    plot(y(:,i)); axis('tight');
    set_graphic_sizes([], 20);
    title(strcat('Micro #',num2str(i)));
end

options.n = n;
w = 128;   % size of the window
q = w/4;    % overlap of the window

clf; X = []; Y = [];
for i=1:s
    X(:,:,i) = perform_stft(x(:,i),w,q, options);
    subplot(s,1,i);
    plot_spectrogram(X(:,:,i));
    set_graphic_sizes([], 20);
    title(strcat('Source #',num2str(i)));
end

exo1()

%% Insert your code here.

mf = size(Y,1);
mt = size(Y,2);
P = reshape(Y, [mt*mf p]);
P = [real(P);imag(P)];

npts = 6000;

sel = randperm(n); sel = sel(1:npts);
clf;
plot(y(sel,1), y(sel,2), '.');
axis([-1 1 -1 1]*5);
set_graphic_sizes([], 20);
title('Time domain');

exo2()

%% Insert your code here.

Theta = mod(atan2(P(:,2),P(:,1)), pi());

nbins = 100;
[h,t] = hist(Theta, nbins);
h = h/sum(h);
clf;
bar(t,h); axis('tight');

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.

A = reshape(Y, [mt*mf p]);

C = abs( M1'*A' );

[tmp,I] = compute_max(C, 1);
I = reshape(I, [mf mt]);

T = .05;
D = sqrt(sum( abs(Y).^2, 3));
I = I .* (D>T);

clf;
imageplot(I(1:mf/2,:));
axis('normal');
set_colormap('jet');

Proj = M1'*A';
Xr = [];
for i=1:s
    Xr(:,:,i) = reshape(Proj(i,:), [mf mt]) .* (I==i);
end

for i=1:s
    xr(:,i) = perform_stft(Xr(:,:,i),w,q, options);
end

clf;
for i=1:s
    subplot(s,1,i);
    plot(xr(:,i)); axis('tight');
    set_graphic_sizes([], 20);
    title(strcat('Estimated source #',num2str(i)));
end

i = 1;
sound(x(:,i),fs);
sound(xr(:,i),fs);
