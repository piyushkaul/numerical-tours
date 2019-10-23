
warning off
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/fastmarching_0_implementing')
warning on

n = 40;

neigh = [[1;0] [-1;0] [0;1] [0;-1]];

boundary = @(x)mod(x-1,n)+1;

ind2sub1 = @(k)[rem(k-1, n)+1; (k - rem(k-1, n) - 1)/n + 1]; 
sub2ind1 = @(u)(u(2)-1)*n+u(1);
Neigh = @(k,i)sub2ind1( boundary(ind2sub1(k)+neigh(:,i)) );

W = ones(n);

x0 = [n/2;n/2];

I = sub2ind1(x0);

D = zeros(n)+Inf; 
D(I) = 0;

S = zeros(n);
S(I) = 1;

[tmp,j] = sort(D(I)); j = j(1);
i = I(j); I(j) = [];

S(i) = -1;

J = [Neigh(i,1); Neigh(i,2); Neigh(i,3); Neigh(i,4)];

J(S(J)==-1) = [];

J1 = J(S(J)==0);
I = [I; J1];
S(J1) = 1;

for j=J'
    dx = min( D([Neigh(j,1) Neigh(j,2)]) );
    dy = min( D([Neigh(j,3) Neigh(j,4)]) );
    D(j) = min(dx+W(j), dy+W(j));
end

exo1()

%% Insert your code here.

displ = @(D)cos(2*pi*5*D/max(D(:)));
clf;
imageplot(displ(D));
colormap jet(256);

D(j) = min(dx+W(j), dy+W(j));

Delta = 2*W(j) - (dx-dy)^2;
if Delta>=0
    D(j) = (dx+dy+sqrt(Delta))/2;
else
    D(j) = min(dx+W(j), dy+W(j));
end

exo2()

%% Insert your code here.

clf;
imageplot(displ(D));
colormap jet(256);

n = 100;
x = linspace(-1,1,n);
[Y,X] = meshgrid(x,x);
sigma = .2;
W = 1 + 8 * exp(-(X.^2+Y.^2)/(2*sigma^2));

clf;
imageplot(W);

x0 = round([.1;.1]*n); %  [.8;.8]]*n);

exo3()

%% Insert your code here.

options.order = 2;
G0 = grad(D, options);

G = G0 ./ repmat( sqrt( sum(G0.^2, 3) ), [1 1 2]);

tau = .8;

x1 = round([.9;.88]*n);
gamma = x1;

Geval = @(G,x)[interp2(1:n,1:n,G(:,:,1),x(2),x(1)); ...
             interp2(1:n,1:n,G(:,:,2),x(2),x(1)) ];

g = Geval(G, gamma(:,end));

gamma(:,end+1) = gamma(:,end) - tau*g;

exo4()

%% Insert your code here.

clf; hold on;
imageplot(W); colormap gray(256);
h = plot(gamma(2,:),gamma(1,:), '.b'); set(h, 'LineWidth', 2);
h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);
axis ij;
