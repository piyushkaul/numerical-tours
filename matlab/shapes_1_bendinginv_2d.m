
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/shapes_1_bendinginv_2d')

clear options;
q = 400;
name = 'centaur1';
S = load_image(name,q);
S = perform_blurring(S,5);
S = double( rescale( S )>.5 );
if S(1)==1 
    S = 1-S;
end

b = compute_shape_boundary(S);
L = size(b,2);

lw = 2;
clf; hold on;
imageplot(-S);
plot(b(2,:), b(1,:), 'r', 'LineWidth', lw); axis ij;

S1 = perform_convolution(S, ones(3)/9)<.01; S1=double(S1==0);
C = zeros(q)-Inf; C(S1==1) = +Inf;
options.constraint_map = C;

geod = @(x)perform_fast_marching(ones(q), x, options);

x = [263; 55];
options.nb_iter_max = Inf;
U = geod(x);

clf;
options.display_levelsets = 1;
options.pstart = x;
options.nbr_levelsets = 60;
U(S==0) = Inf;
display_shape_function(U, options);

exo1()

%% Insert your code here.

N = 1000;

N0 = round(.4*N);

I = round(linspace(1,L+1,N0+1));
X = round(b(:,I(1:end-1)));

[y,x] = meshgrid(1:q,1:q);
I = find(S==1);
I = I(randperm(length(I))); I = I(1:N-N0);
X(:,end+1:N) = [x(I),y(I)]';

clf; hold on;
imageplot(1-S);
plot(X(2,:), X(1,:), 'r.', 'MarkerSize', 15);
axis('ij');

exo2()

%% Insert your code here.

clf;
imageplot(delta(1:N0,1:N0));
colormap(jet(256));

J = eye(N) - ones(N)/N;

K = -1/2 * J*(delta.^2)*J;

opt.disp = 0; 
[Y, v] = eigs(K, 2, 'LR', opt);
Y = Y .* repmat(sqrt(diag(v))', [N 1]);
Y = Y';

theta = -.8*pi;
uv = Y(:,1:N0);
uv = [cos(theta)*uv(1,:) + sin(theta)*uv(2,:); - sin(theta)*uv(1,:) + cos(theta)*uv(2,:)];

clf;
h = plot(uv(2,:), uv(1,:)); 
axis('ij'); axis('equal'); axis('off');
set(h, 'LineWidth', lw);

Stress = @(d)sqrt( sum( abs(delta(:)-d(:)).^2 ) / N^2 );

D = @(Y)sqrt( repmat(sum(Y.^2),N,1) + repmat(sum(Y.^2),N,1)' - 2*Y'*Y);

Y = X/q;

remove_diag = @(b)b - diag(sum(b));
B = @(D1)remove_diag( -delta./max(D1,1e-10) );

Y = Y * B(D(Y))' / N;

exo3()

%% Insert your code here.

clf;
plot(s, '.-', 'LineWidth', 2, 'MarkerSize', 20);
axis('tight');

clf;
h = plot(Y(2,[1:N0 1]), Y(1,[1:N0 1])); 
axis('ij'); axis('equal'); axis('off');
set(h, 'LineWidth', lw);

exo4()

%% Insert your code here.
