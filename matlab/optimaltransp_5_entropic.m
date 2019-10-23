
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/optimaltransp_5_entropic')

N = [300,200];

d = 2;

x = rand(2,N(1))-.5;

theta = 2*pi*rand(1,N(2));
r = .8 + .2*rand(1,N(2));
y = [cos(theta).*r; sin(theta).*r];

plotp = @(x,col)plot(x(1,:)', x(2,:)', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', 2);

clf; hold on;
plotp(x, 'b');
plotp(y, 'r');
axis('off'); axis('equal');

x2 = sum(x.^2,1); y2 = sum(y.^2,1);
C = repmat(y2,N(1),1)+repmat(x2.',1,N(2))-2*x.'*y;

p = ones(N(1),1)/N(1);
q = ones(N(2),1)/N(2);

gamma = .01;

xi = exp(-C/gamma);

b = ones(N(2),1);

a = p ./ (xi*b);
b = q ./ (xi'*a);

exo1()

%% Insert your code here.

Pi = diag(a)*xi*diag(b);

clf;
imageplot(Pi);

exo2()

%% Insert your code here.

Pi = diag(a)*xi*diag(b);

clf;
hold on;
A = sparse( Pi .* (Pi> min(1./N)*.7) ); [i,j,~] = find(A);
h = plot([x(1,i);y(1,j)], [x(2,i);y(2,j)], 'k');
set(h, 'LineWidth', 2); % weaker connections.
A = sparse( Pi .* (Pi> min(1./N)*.3) ); [i,j,~] = find(A);
h = plot([x(1,i);y(1,j)], [x(2,i);y(2,j)], 'k:');
set(h, 'LineWidth', 1);
plotp(x, 'b'); % plot the two point clouds.
plotp(y, 'r');
axis('off'); axis('equal');

N = 200;

t = (0:N-1)'/N;

Gaussian = @(t0,sigma)exp( -(t-t0).^2/(2*sigma^2) );
normalize = @(p)p/sum(p(:));
sigma = .06;
p = Gaussian(.25,sigma); 
q = Gaussian(.8,sigma);

vmin = .02;
p = normalize( p+max(p)*vmin);
q = normalize( q+max(q)*vmin);

clf;
subplot(2,1,1);
bar(t, p, 'k'); axis tight;
subplot(2,1,2);
bar(t, q, 'k'); axis tight;

gamma = (.03)^2;

[Y,X] = meshgrid(t,t);
xi = exp( -(X-Y).^2 / gamma);

b = ones(N,1);

a = p ./ (xi*b);
b = q ./ (xi'*a);

exo3()

%% Insert your code here.

Pi = diag(a)*xi*diag(b);
clf;
imageplot(log(Pi+1e-5));

s = (xi*(b.*t)) .* a ./ p;

clf; hold on;
imagesc(t,t,log(Pi+1e-5)); colormap gray(256);
plot(s,t, 'r', 'LineWidth', 3);
axis image; axis off; axis ij;

N = 70;

names = {'disk' 'twodisks' 'letter-x' 'letter-z'};
r = [.35 .19 .12*N .12*N]; 
vmin = .05;
P = [];
for i=1:length(names)
    options.radius = r(i);
    p = load_image(names{i},N, options); 
    p = normalize( rescale(p)+vmin );
    P(:,:,i) = p;
end
K = size(P,3);

a = mat2cell(P, N,N,ones(K,1));
clf;
imageplot(a, '', 2,2);

gamma = (.04)^2;

n = 41; % width of the convolution kernel
t = linspace(-n/(2*N),n/(2*N),n)';
g = exp(-t.^2 / gamma); g2 = g*g';
xi = @(x)conv2(conv2(x, g, 'same')', g, 'same')';

clf;
imageplot({P(:,:,1) xi(P(:,:,1))});

lambda = ones(K,1)/K;

b = ones(N,N,K); a = b;

for k=1:K
    a(:,:,k) = P(:,:,k) ./ xi(b(:,:,k));
end

q = zeros(N);
for k=1:K
    q = q + lambda(k) * log( max(1e-19, b(:,:,k) .* xi(a(:,:,k)) ) );
end
q = exp(q);

for k=1:K
    b(:,:,k) = q ./ xi(a(:,:,k));
end

exo4()

%% Insert your code here.

clf;
imageplot(q);

exo5()

%% Insert your code here.
