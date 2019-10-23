
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/meshproc_5_pde')

name = 'bunny';
clear options;
options.name = name;
[X0,F] = read_mesh(name);
n = size(X0,2);
options.name = name;

clf;
plot_mesh(X0,F,options);

exo1()

%% Insert your code here.

d = full( sum(W,1) );
D = spdiags(d(:), 0, n,n);
L = D - W;

iD = spdiags(d(:).^(-1), 0, n,n);
tW = iD * W;
tL = iD * L;

M = load_image('lena',256);

v = X0 - repmat(mean(X0,2), [1 n]);
theta = acos(v(1,:)./sqrt(sum(v.^2)))/pi;
phi = (atan2(v(2,:),v(3,:))/pi+1)/2;

x = linspace(0,1,size(M,1));
f0 = interp2(x,x,M',theta,phi)';

options.face_vertex_color = f0(:);
clf;
plot_mesh(X0,F, options);
lighting none;

Tmax = 50;

tau = .4;

niter = ceil(Tmax/tau);

f = f0;

f = f - tau*tL*f;

exo2()

%% Insert your code here.

X = X0;

X = X - tau*X*tL';

exo3()

%% Insert your code here.

n = size(X0,2);

sigma = .002;

npoints = 16;

f0 = zeros(n,1);

D = zeros(n,1) + 1e10;
for i=1:npoints
    % select a point farthest away
    [tmp,p] = max(D);
    % compute the distance beetween the point and the other vertices
    d = sum( (X0 - X0(:,p)*ones(1,n)).^2 );
    f0 = f0 + (-1)^mod(i,2) *exp( -d(:)/(2*sigma.^2) );
    D = min(D,d(:));
end

options.face_vertex_color = f0;
clf;
plot_mesh(X0,F, options);
colormap(jet(256));

Tmax = 80;
tau = .5;
niter = round(Tmax/tau);

f1 = f0;

update = @(f,f1)deal(2*f - f1 - tau^2 * tL*f, f);
[f,f1] = update(f,f1);

exo4()

%% Insert your code here.
