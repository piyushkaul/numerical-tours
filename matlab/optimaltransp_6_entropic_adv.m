
addpath('solutions/optimaltransp_6_entropic_adv')

n = 100;
m = 200; 
d = 2; % Dimension of the clouds.
a = ones(n,1)/n; 
b = ones(1,m)/m;

x = rand(2,n)-.5;
theta = 2*pi*rand(1,m);
r = .8 + .2*rand(1,m);
y = [cos(theta).*r; sin(theta).*r];

plotp = @(x,col)plot(x(1,:)', x(2,:)', 'o', 'MarkerSize', 9, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', 2);
clf; hold on;
plotp(x, 'b');
plotp(y, 'r');
axis('off'); axis('equal');

distmat = @(x,y)sum(x.^2,1)' + sum(y.^2,1) - 2*x.'*y;
C = distmat(x,y);

mina = @(H,epsilon)-epsilon*log( sum(a .* exp(-H/epsilon),1) );
minb = @(H,epsilon)-epsilon*log( sum(b .* exp(-H/epsilon),2) );

minaa = @(H,epsilon)mina(H-min(H,[],1),epsilon) + min(H,[],1);
minbb = @(H,epsilon)minb(H-min(H,[],2),epsilon) + min(H,[],2);

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

z = x; % initialization

epsilon = .01;
niter = 300;
C = distmat(z,y); 
for it=1:niter
    g = mina(C-f,epsilon);
    f = minb(C-g,epsilon);
end
P = a .* exp((f+g-C)/epsilon) .* b;

G = z - ( y*P' ) ./ a';

clf; hold on;
plotp(x, 'b');
plotp(y, 'r');
for i=1:n
    plot([x(1,i), x(1,i)-G(1,i)], [x(2,i), x(2,i)-G(2,i)], 'k');
end
axis('off'); axis('equal');

tau = .1;

z = z - tau * G;

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.

z = randn(2,n)*.2;
y = randn(2,m)*.2; y(1,:) = y(1,:)*.05 + 1;

A = eye(2); h = zeros(2,1);

clf; hold on;
plotp(A*z+h, 'b'); plotp(y, 'r'); 
axis('off'); axis('equal');

x = A*z+h;
C = distmat(x,y); 
f = zeros(n,1);
for it=1:niter
	g = mina(C-f,epsilon);
	f = minb(C-g,epsilon);
end
P = a .* exp((f+g-C)/epsilon) .* b;

v = a' .* z - ( y*P' );

nabla_A = v*z';
nabla_h = sum(v,2);

exo5()

%% Insert your code here.
