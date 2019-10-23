
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/inverse_4_inpainting_variational')

name = 'cameraman';
n = 256;
f0 = rescale( load_image(name, n) );

clf;
imageplot(f0);

rho = .7;

Gamma = rand(n)>rho;

Phi = @(f)f.*Gamma;

y = Phi(f0);

clf;
imageplot(y);

Pi = @(f)f.*(1-Gamma) + y.*Gamma;

Delta = @(f)div(grad(f));

tau = .8/4;

exo1()

%% Insert your code here.

clf;
plot(E); axis('tight');
set_label('Iteration #', 'E');

clf;
imageplot(f, strcat(['Inpainted, SNR=' num2str(snr(f0,f),3) 'dB']));

epsilon = 1e-2;

Amplitude = @(u)sqrt(sum(u.^2,3)+epsilon^2);
Neps = @(u)u./repmat(Amplitude(u), [1 1 2]);

tau = .9*epsilon/4;

G = @(f)-div(Neps(grad(f)));

exo2()

%% Insert your code here.

clf;
imageplot(clamp(f), strcat(['SNR=' num2str(snr(f0,f),3) 'dB']));

clf;
plot(J); 
axis('tight');
set_label('Iteration #', 'J_\epsilon');

n = 256;
f0 = load_image('parrot', n);
f0 = rescale( sum(f0,3) );

clf;
imageplot(f0);

Gamma = load_image('parrot-mask', n);
Gamma = double(rescale(Gamma)>.5);

Phi = @(f)f.*Gamma;

y = Phi(f0);

clf;
imageplot(y);

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.
