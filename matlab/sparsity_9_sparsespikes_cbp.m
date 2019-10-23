
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/sparsity_9_sparsespikes_cbp')

sigma = .03;
phi  = @(t)exp(-t.^2/(2*sigma^2));
phi1 = @(t)-t/(sigma^2).*exp(-t.^2/(2*sigma^2));

N = 32; 
Delta = 1/N;

rho = 64;
P = N*rho;

convol  = @(x,h)real(ifft(fft(x).*fft(h)));
convolS = @(x,h)real(ifft(fft(x).*conj(fft(h))));

upsample = @(x,rho)upsampling(x,1,rho);
downsample = @(y,rho)y(1:rho:end);

t = [0:P/2, -P/2+1:-1]' / P;
t1 = (0:P-1)'/P;

phi_d = phi(t);
phi1_d = phi1(t);

Phi  = @(x)convol(upsample(x,rho),phi_d);
PhiS = @(x)downsample(convolS(x,phi_d),rho);
Psi  = @(s)convol(upsample(s,rho),phi1_d)*Delta/2;
PsiS = @(s)downsample(convolS(s,phi1_d),rho)*Delta/2;
Gamma  = @(u)Phi(u(:,1)) - Psi(u(:,2));
GammaS = @(y)[PhiS(y), -PsiS(y)];

k = 3; % number of spikes
kappa = .9;
I = round( N/(2*k):N/k:N );
a0 = zeros(N,1); a0(I) = [.6 1 .8];
d0 = zeros(N,1); d0(I) = [-.2 1 -.7] * kappa;
b0 = d0.*a0;
x0 = (0:N-1)'/N + d0*Delta/2;

y = zeros(P,1);
for i=1:length(x0)
    T = t-x0(i); T = mod(T,1); T(T>.5) = T(T>.5)-1;
    y = y + a0(i) * phi( T );
end

y0 = Phi(a0);

y1 = Gamma([a0 b0]);

lw = 2; msB = 30;
mystem = @(x,y, col)stem(x, y, [col '.--'], 'MarkerSize', msB, 'LineWidth', lw);
subplot(3,1,1); hold on;
mystem(x0(I), a0(I), 'k'); box on;
plot(t1, y, 'LineWidth', lw); axis tight; title('Observations');
subplot(3,1,2);
plot(t1, [y-y0 y-y1], 'LineWidth', lw); axis tight; title('0th order');
legend('0th order', '1st order');
subplot(3,1,3);
plot(t1, y-y1, 'g', 'LineWidth', lw); axis tight; title('1st order');

C  = @(u)u(:,2)+1i*u(:,1);
Ci = @(v)[imag(v), real(v)];
ProjOct = @(v)max(real(v),0) + 1i*max(imag(v),0);
ProjC = @(u)Ci(exp(1i*pi/4)*ProjOct(exp(-1i*pi/4)*C(u)));

ProxJ = @(w,lambda)ProjC( w-[lambda*ones(size(w,1),1) zeros(size(w,1),1)] );

exo1()

%% Insert your code here.

lambda = 40;

F = @(u)1/2*norm(y-Gamma(u), 'fro')^2;
gradF = @(u)GammaS(Gamma(u)-y);

E = @(u)F(u) + lambda*norm(u(:,1),1);

u = zeros(N,2);

exo2()

%% Insert your code here.

a = u(:,1); 
b = u(:,2);
delta = Delta/2 * b./a; delta(a<1e-9) = 0;
x = (0:N-1)'/N + delta;

J = find(a>1e-3);
t = (0:N-1)'/N;
s = (0:P-1)'/P;
clf; hold on;
plot(s, y, 'LineWidth', 2);
mystem(x0(I), a0(I), 'k'); % initial spikes 
mystem(t(J) + delta(J), a(J), 'r');  % recovered spikes
axis([0 1 0 1]);
box on;

subplot(2,1,1);
hold on;
mystem(1:N, a0, 'k');
plot(a, 'r.', 'MarkerSize', 20); 
axis tight; box on; title('a');
subplot(2,1,2);
hold on;
mystem(1:N, b0, 'k');
plot(b, 'r.', 'MarkerSize', 20); 
axis tight; box on; title('b');

exo3()

%% Insert your code here.
