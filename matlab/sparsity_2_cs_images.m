
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/sparsity_2_cs_images')

name = 'boat';
n = 256;
f = load_image(name, n);
f = rescale(f);

k0 = 2;
J = log2(n)-k0;
Wav  = @(f)perform_wavelet_transf(f,J,+1);
WavI = @(x)perform_wavelet_transf(x,J,-1);

fw = Wav(f);

clf;
plot_wavelet(fw, J);

exo1()

%% Insert your code here.

A = ones(n,n); A(1:2^J,1:2^J) = 0;
I0 = find(A==1);
x0 = fw(I0);

N = length(x0);

P0 = (n/2^k0)^2;

P = 4 * P0;

sigma1 = randperm(N)'; 
sigma2 = randperm(N)'; 
S1 = @(x)x(sigma1);
S2 = @(x)x(sigma2);

sigma1S = 1:N; sigma1S(sigma1) = 1:N; 
sigma2S = 1:N; sigma2S(sigma2) = 1:N; 
S1S = @(x)x(sigma1S);
S2S = @(x)x(sigma2S);

downarrow = @(x)x(1:P); 
Phi = @(x)downarrow(S2(dct(S1(x))));

uparrow = @(x)[x; zeros(N-P,1)];
PhiS = @(x)S1S(idct(S2S(uparrow(x))));

y = Phi(x0);

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

mu = 1;
gamma = 1;

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.

w = 4;

v = 1:w:n;
dv = 0:w-1;
[dX,dY,X,Y] = ndgrid(dv,dv,v,v);
q = size(X,3);
dX = reshape(dX, [w*w q*q]);
dY = reshape(dY, [w*w q*q]);
X = reshape(X, [w*w q*q]);
Y = reshape(Y, [w*w q*q]);

I = find( sum(X+dX>n | Y+dY>n)  );
X(:,I) = [];
Y(:,I) = [];
dX(:,I) = [];
dY(:,I) = [];

U = zeros(n,n);
U(I0) = 1:N;
Ind = X+dX + (Y+dY-1)*n;
I = U(Ind);

I(:,sum(I==0)>0) = [];

G = @(x)sum( sqrt(sum(x(I).^2)) );

[A,tmp] = meshgrid( randperm(size(I,2)) , ones(w*w,1));
x = zeros(N,1); x(I) = A;
Z = zeros(n,n); Z(I0) = x;
clf; 
imageplot(Z); 
colormap jet(256);

exo6()

%% Insert your code here.

exo7()

%% Insert your code here.

exo8()

%% Insert your code here.
