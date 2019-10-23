
warning off
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/denoisingsimp_4_denoiseregul')
warning on

n = 256;
name = 'hibiscus';
f0 = load_image(name,n);
f0 = rescale( sum(f0,3) );

Gr = grad(f0);

d = sqrt(sum3(Gr.^2,3));

clf;
imageplot(Gr, strcat(['grad']), 1,2,1);
imageplot(d, strcat(['|grad|']), 1,2,2);

sob = sum(d(:).^2);

sigma = .1;
y = f0 + randn(n,n)*sigma;

yF = fft2(y);

x = [0:n/2-1, -n/2:-1];
[Y,X] = meshgrid(x,x);
S = (X.^2 + Y.^2)*(2/n)^2;

lambda = 20;

fSob = real( ifft2( yF ./ ( 1 + lambda*S) ) );

clf;
imageplot(clamp(fSob));

exo1()

%% Insert your code here.

esob = snr(f0,fSob0);  enoisy = snr(f0,y);
clf;
imageplot(clamp(y), strcat(['Noisy ' num2str(enoisy,3) 'dB']), 1,2,1);
imageplot(clamp(fSob0), strcat(['Sobolev regularization ' num2str(esob,3) 'dB']), 1,2,2);

exo2()

%% Insert your code here.

u = linspace(-5,5)';
clf;
subplot(2,1,1); hold('on');
plot(u, abs(u), 'b');
plot(u, sqrt(.5^2+u.^2), 'r');
title('\epsilon=1/2'); axis('square');
subplot(2,1,2); hold('on');
plot(u, abs(u), 'b');
plot(u, sqrt(1^2+u.^2), 'r');
title('\epsilon=1'); axis('square');

epsilon_list = [1e-9 1e-2 1e-1 .5];
clf;
for i=1:length(epsilon_list)
    G = div( Gr ./ repmat( sqrt( epsilon_list(i)^2 + d.^2 ) , [1 1 2]) );
    imageplot(G, strcat(['epsilon=' num2str(epsilon_list(i))]), 2,2,i);
end

epsilon = 1e-2;

lambda = .1;

tau = 2 / ( 1 + lambda * 8 / epsilon);

fTV = y;

Gr = grad(fTV);
d = sqrt(sum3(Gr.^2,3));
G = -div( Gr ./ repmat( sqrt( epsilon^2 + d.^2 ) , [1 1 2]) );

fTV = fTV - tau*( y-fTV + lambda* G);

exo3()

%% Insert your code here.

clf;
imageplot(fTV);

exo4()

%% Insert your code here.

etvr = snr(f0,fTV0); 
clf;
imageplot(clamp(y), strcat(['Noisy ' num2str(enoisy,3) 'dB']), 1,2,1);
imageplot(clamp(fTV0), strcat(['TV regularization ' num2str(etvr,3) 'dB']), 1,2,2);

exo5()

%% Insert your code here.
