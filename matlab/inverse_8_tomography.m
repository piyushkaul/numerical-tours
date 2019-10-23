
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/inverse_8_tomography')

n = 128;

name = 'phantom';
f0 = load_image(name, n);

t = linspace(-1,1,n);
[Y,X] = meshgrid(t,t);
rotation = @(f,t)interp2(Y,X,f, sin(t)*X + cos(t)*Y, cos(t)*X - sin(t)*Y, 'cubic', 0);

clf;
imageplot(abs(rotation(rotation(f0,.2),-.2)-f0));

t = [0:n/2-1 0 -n/2+1:-1]';
[Y,X] = meshgrid(t,t);
shearx = @(f,a)real( ifft( fft(f) .* exp(-a*2i*pi*X.*Y/n) ) );
shearx = @(f,a)fftshift(shearx(fftshift(f),a));
sheary = @(f,a)shearx(f',a)';

exo1()

%% Insert your code here.

rotation = @(f,t)shearx( sheary( shearx(f,-tan(t/2)) ,sin(t)) ,-tan(t/2));

exo2()

%% Insert your code here.

rotation = @(f,t)rotation(rotation(f,t/2),t/2);
rotation = @(f,t)rotation(rotation(f,t/2),t/2);

exo3()

%% Insert your code here.

theta = .2;
e = norm( rotation(rotation(f0,theta),-theta)-f0 );
disp(['Error (should be 0): ' num2str(e, 2) '.']);

theta = .1;
a = randn(n); b = randn(n);
dotp = @(a,b)sum(a(:).*b(:));
e = dotp(a,rotation(b,theta)) - dotp(b,rotation(a,-theta));
disp(['Error (should be 0): ' num2str(e, 2) '.']);

m = 64;

tlist = linspace(-pi/2,pi/2, m+1); tlist(end) = [];

Phi  = @(f)perform_radon_transform(f, tlist, +1, rotation);
PhiS = @(R)perform_radon_transform(R, tlist, -1, rotation);

clf;
imageplot(Phi(f0));

a = randn(n,m); b = randn(n);
e = dotp(a,Phi(b)) - dotp(PhiS(a),b);
disp(['Error (should be 0): ' num2str(e, 2) '.']);

m = 16;
tlist = linspace(-pi/2,pi/2, m+1); tlist(end) = [];
Phi  = @(f)perform_radon_transform(f, tlist, +1, rotation);
PhiS = @(R)perform_radon_transform(R, tlist, -1, rotation);

y = Phi(f0);

clf;
plot(y); axis tight;

exo4()

%% Insert your code here.

return;

Amplitude = @(u)sqrt(sum(u.^2,3));
F1 = @(u)sum(sum(Amplitude(u)));
u = @(z)reshape(z(1:n*m),n,m);
v = @(z)reshape(z(n*m+1:end),n,n,2);

ProxF1 = @(u,lambda)max(0,1-lambda./repmat(Amplitude(u), [1 1 2])).*u;
ProxF = @(z,lambda)[y(:); reshape(ProxF1(v(z),lambda),n*n*2,1) ];
ProxFS = @(y,sigma)y-sigma*ProxF(y/sigma,1/sigma);
ProxG = @(x,lambda)x;
K  = @(f)[reshape(Phi(f),n*m,1); reshape(grad(f), n*n*2,1)];
KS = @(z)PhiS(u(z)) - div(v(z));

L = n*m;
sigma = 10;
tau = .9/(L*sigma);
theta = 1;

f = fL2;

g = K(f)*0;
f1 = f;

niter = 100;
E = []; C = []; F = [];
for i=1:niter    
    % update
    fold = f;
    g = ProxFS( g+sigma*K(f1), sigma);
    f = ProxG(  f-tau*KS(g), tau);
    f1 = f + theta * (f-fold);
    % monitor the decay of the energy
    E(i) = F1(grad(f));
    F(i) = norm(y-Phi(f), 'fro');
    C(i) = snr(f0,f);
end
clf;
h = plot(E);
set(h, 'LineWidth', 2);
axis('tight');
