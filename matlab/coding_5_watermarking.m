
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/coding_5_watermarking')

n = 256;
name = 'hibiscus';
f = load_image(name, n);
f = rescale(sum(f,3));

clf;
imageplot(f);

Jmin = log2(n)-2;
Psi  = @(f)perform_wavelet_transf(f, Jmin, +1);
PsiS = @(a)perform_wavelet_transf(a, Jmin, -1);

a = Psi(f);

clf;
plot_wavelet(a,Jmin);

A = ones(n); A(1:2^Jmin,1:2^Jmin) = 0;
I = find(A(:));
P = length(I);

x0 = a(I);

w = randn(P,1);

psnr_embedding = 50;

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

x = x0 + rho*abs(x0).*w;

disp(['PSNR(x,x0) = ' num2str(psnr(x,x0,1), 3) 'dB.']);

a1 = a; a1(I) = x;
f1 = PsiS(a1);

delta = f-f1;
clf;
imageplot( clamp(delta/std(delta(:)),-3,3) );

C = @(y,w)sum(w.*y)./sqrt( sum(w.^2).*sum(y.^2) );

exo3()

%% Insert your code here.

pfa = 1e-3;
T = sqrt(2)/2 * sigma0 * erfinv(1-2*pfa);

exo4()

%% Insert your code here.

tau = .2;

Quant = @(x)floor(abs(x/tau)).*sign(x);
DeQuant = @(v)sign(v) .* (abs(v)+.5) * tau;
A = @(x)DeQuant(Quant(x));

t = linspace(-2,2,500);
plot(t, A(t));
axis('equal');

y = A(x);

a1 = a; a1(I) = y;
f1 = PsiS(a1);

clf;
imageplot(clamp(f1));

disp(['C(y,w) = ' num2str(C(y,w), 2) '.']);

exo5()

%% Insert your code here.

exo6()

%% Insert your code here.

exo7()

%% Insert your code here.

exo8()

%% Insert your code here.
