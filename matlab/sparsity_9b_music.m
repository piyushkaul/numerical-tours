
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/sparsity_9b_music')

ms = 20; % markersize

M = 18;

Q = 1024*4;
z = (0:Q)'/Q;

x0 = [.1 .4 .5 .7 .9]';

N = 5;

a0 = [.7 -.8 .9 1 -.9]';

Phi = @(x)exp(-2i*pi*(0:M)'*x(:)');

sigma = .15;
w = (randn(M+1,1)+1i*randn(M+1,1)); w = w/norm(w);
y0 = Phi(x0)*a0; 
y = y0 + sigma*norm(y0)*w;

PhiT = @(x)exp(2i*pi*x(:)*(-M:M)) / (2*M+1);
f = real( PhiT(z) * [conj(y(end:-1:2)); y]);

clf; hold on;
plot(z, f);
stem(x0,a0, 'r');

L = M/2;

MusicHankel = @(y)hankel(y(1:L),y(L:M));

[U,S,V] = svd(MusicHankel(y0),0); 
S = diag(S);

clf;
plot(S, '.-', 'MarkerSize', ms);
axis tight;

Ubot = U(:,N+1:end);
d = Ubot'*exp(-2i*pi*(0:L-1)'*z(:)');
d = sum(abs(d).^2) / L;

clf; hold on;
plot(z, d, 'b');
stem(x0, 1+x0*0, 'r.');
axis([0 1 0 1]); box on;

B = [];
for j=1:size(Ubot,2)
    u = Ubot(:,j);
    v = flipud(conj(u)); 
    B(:,j) = conv(u,v);
end
C = sum(B,2);

R = roots(C(end:-1:1));

R1 = R(abs( abs(R)-1 )<1e-4);

clf; hold on;
plot(exp(2i*pi*z), 'k');
plot(R, 'b.', 'MarkerSize', ms);
plot(R1, 'r.', 'MarkerSize', ms);
axis equal; box on;
axis([-1 1 -1 1]*2);

[x1,I] = sort(mod(angle(R1),2*pi)/(2*pi));
x1 = x1(1:2:end);

disp( ['Error (should be 0): ' num2str(norm(x1-x0))] );

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.
