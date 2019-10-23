
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/graphics_1_synthesis_gaussian')

n = 512;
name = 'wood';
f = rescale( load_image(name, n) );

clf;
imageplot(f);

z = zeros(1,1,size(f,3));
fe = [z, f(1,:,:), z; f(:,1,:), f, f(:,end,:); z, f(end,:,:), z];

laplacian = @(x)4*x - ( circshift(x,[0 1]) + circshift(x,[1 0]) + circshift(x,[-1 0]) + circshift(x,[0 -1]) );
d = laplacian(fe);
d = d(2:end-1, 2:end-1, :);

[X Y] = meshgrid(0:n-1, 0:n-1);
U = 4 - 2*cos(2.*X*pi/n) - 2*cos(2.*Y*pi/n);

P = fft2(d)./repmat(U, [1 1 size(f,3)]);
P(1,1,:) = sum(sum(f,1),2);
p = real(ifft2(P));

mydisp = @(x)[x x; x x];
clf;
imageplot(mydisp(f), 'Original, periodized',1,2,1);
imageplot(mydisp(p), 'Periodic layer, periodized',1,2,2);

exo1()

%% Insert your code here.

f0 = mean(p,3);

u = f0; u(1,1,:)=0; u(2,1,:)=1;
clf;
imageplot(clamp(u));

w = randn(n)/n; 
w = w-mean(w(:))+1/n^2;

f = real(ifft2(fft2(f0).*fft2(w)));

clf;
u = f0; u(1,1,:)=0; u(2,1,:)=1;
imageplot(clamp(u), 'Input', 1,2,1);
u = f; u(1,1,:)=0; u(2,1,:)=1;
imageplot(clamp(u), 'Synthesized', 1,2,2);

exo2()

%% Insert your code here.

f0 = p;

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.

n0 = n*2;

a = 10/n;

phi = @(t)sin(t/a*pi/2).^2;
psi = @(t)phi(t).*(t<a) + (t>=a & t<=1-a) + phi(1-t).*(t>1-a);

t = linspace(0,1,n)';
g = repmat( psi(t)*psi(t)', [1 1 size(f0,3)]);

mu0 = repmat(mean(mean(f0,1),2), [n n 1]);
f1 = repmat(mean(mean(f0,1),2), [n0 n0 1]);
f1(end/2-n/2+1:end/2+n/2, end/2-n/2+1:end/2+n/2, :) = mu0 + n0/n * g .* (f0-mu0);

clf;
u = f1; u(1,1,:)=0; u(2,1,:)=1;
imageplot(clamp(u));

exo5()

%% Insert your code here.
