
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/sparsity_8_sparsespikes_measures')

ms = 20;
lw = 1;

P = 2048*8;
options.P = P;
u = (0:P-1)'/P;

fc = 6;
N = 2*fc+1;

Fourier = @(fc,x)exp(-2i*pi*(-fc:fc)'*x(:)');

Phi  = @(fc,x,a)Fourier(fc,x)*a;
PhiS = @(fc,u,p)real( Fourier(fc,u)'* p );

delta = .7/fc;

x0 = [.5-delta .5 .5+delta]';
a0 = [1 1 -1]';
n = length(x0);

y0 = Phi(fc,x0,a0);

sigma = .12 * norm(y0);
w = fftshift( fft(randn(N,1)) ); 
w = w/norm(w) * sigma;
y = y0 + w;

f0 = PhiS(fc,u,y0);
f = PhiS(fc,u,y);
clf; hold on;
plot(u, [f0 f]);
stem(x0, 10*sign(a0), 'k.--', 'MarkerSize', ms, 'LineWidth', 1);
axis tight; box on;
legend('\Phi^* y_0', '\Phi^* y');

lambda = 1;

dotp = @(x,y)real(x'*y);
Xmat = @(p,Q)[Q, p; p', 1];
Qmat = @(X)X(1:end-1,1:end-1);
pVec = @(X)X(1:end-1,end);

f = @(p)1/2*norm( y/lambda-p )^2;

Proxf = @(p,gamma)( p + gamma*y/lambda )/(1+gamma);

ProxG = @(X,gamma)perform_sdp_projection(X);

ProxF = @(X,gamma)Xmat( Proxf(pVec(X),gamma/2), perform_sos_projection(Qmat(X)) );

rProxF = @(x,tau)2*ProxF(x,tau)-x;
rProxG = @(x,tau)2*ProxG(x,tau)-x;

X = zeros(2*fc+2);

gamma = 1/10;
mu = 1;

Y = X;
ObjVal = []; ConstrSDP = []; ConstrSOS = [];
niter = 300;
for i=1:niter 
    % record energies
    ObjVal(i) = f(pVec(X));
    ConstrSDP(i) = min(real(eig(X)));
    ConstrSOS(i) = norm(perform_sos_projection(Qmat(X))-Qmat(X), 'fro'); 
    % iterate
	Y = (1-mu/2)*Y + mu/2*rProxF( rProxG(Y,gamma),gamma );        
	X = ProxG(Y,gamma);
end
p = pVec(X);

etaLambda = PhiS(fc,u,p);
clf; hold on;
stem(x0, sign(a0), 'k.--', 'MarkerSize', ms, 'LineWidth', lw);
plot([0 1],  [1 1], 'k--', 'LineWidth', lw); 
plot([0 1], -[1 1], 'k--', 'LineWidth', lw);
plot(u, etaLambda, 'b', 'LineWidth', lw);
axis([0 1 -1.1 1.1]);
set(gca, 'XTick', [], 'YTick', [0 1]); box on;

clf; hold on;
stem(x0, sign(abs(a0)), 'k.--', 'MarkerSize', ms, 'LineWidth', lw);
plot([0 1],  [1 1], 'k--', 'LineWidth', lw);
plot(u, 1-abs(etaLambda).^2, 'b', 'LineWidth', lw);
axis([0 1 -.1 1.1]);
set(gca, 'XTick', [], 'YTick', [0 1]);
box on;

c = -conv(p,flipud(conj(p)));
c(N)=1+c(N);

R = roots(flipud(c));

clf;
plot(real(R),imag(R),'*');
hold on;
plot(cos(2*pi*x0), sin(2*pi*x0),'ro');
plot( exp(1i*linspace(0,2*pi,200)), '--' );
hold off;
legend('Roots','Support of x'); 
axis equal; axis([-1 1 -1 1]*1.5);

tol = 1e-2;
R0 = R(abs(1-abs(R)) < tol);

[~,I]=sort(angle(R0));
R0 = R0(I); R0 = R0(1:2:end);

x = angle(R0)/(2*pi);
x = sort(mod(x,1));

Phix = Fourier(fc,x);
s = sign(real(Phix'*p));

a = real(Phix\y - lambda*pinv(Phix'*Phix)*s );

clf; hold on;
stem(x0, a0, 'k.--', 'MarkerSize', ms, 'LineWidth', 1);
stem(x, a, 'r.--', 'MarkerSize', ms, 'LineWidth', 1);
axis([0 1 -1.1 1.1]); box on;
legend('Original', 'Recovered');

exo1()

%% Insert your code here.

d = 1; % number of derivatives
w = ones(N,1);
Gamma = [];
for i=0:d
    Gamma = [Gamma, diag(w) * Fourier(fc,x0)];
    % derivate the filter
    w = w .* 2i*pi .* (-fc:fc)';
end

pV = pinv(Gamma') * [sign(a0); zeros(n,1)];
etaV = PhiS(fc, u, pV);

clf; hold on;
stem(x0, sign(a0), 'k.--', 'MarkerSize', ms, 'LineWidth', lw);
plot([0 1],  [1 1], 'k--', 'LineWidth', lw); 
plot([0 1], -[1 1], 'k--', 'LineWidth', lw);
plot(u, etaV, 'b', 'LineWidth', lw);
axis([0 1 -1.1 1.1]);
set(gca, 'XTick', [], 'YTick', [-1 1]);
box on;

exo2()

%% Insert your code here.
