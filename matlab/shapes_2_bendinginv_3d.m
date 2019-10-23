
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/shapes_2_bendinginv_3d')

name = 'camel';
options.name = name;
[V,F] = read_mesh(name);
N = size(V,2);

clf;
plot_mesh(V,F, options);

i = 1;
U = perform_fast_marching_mesh(V, F, i);

options.method = 'continuous';
J = randperm(N); J = J(1:50);
paths = compute_geodesic_mesh(U, V, F, J, options);

clf;
plot_fast_marching_mesh(V, F, U, paths, options);

exo1()

%% Insert your code here.

d = 3;

J = eye(N) - ones(N)/N;

K = -1/2 * J*(delta.^2)*J;

opt.disp = 0; 
[Y, v] = eigs(K, d, 'LR', opt);
Y = Y .* repmat(sqrt(diag(v))', [N 1]);
Y = Y';

clf;
plot_mesh(Y,F, options);

Stress = @(d)sqrt( sum( abs(delta(:)-d(:)).^2 ) / N^2 );

D = @(Y)sqrt( repmat(sum(Y.^2),N,1) + repmat(sum(Y.^2),N,1)' - 2*Y'*Y);

Y = V;

remove_diag = @(b)b - diag(sum(b));
B = @(D1)remove_diag( -delta./max(D1,1e-10) );

Y = Y * B(D(Y))' / N;

exo2()

%% Insert your code here.

clf;
plot(s, '.-', 'LineWidth', 2, 'MarkerSize', 20);
axis('tight');

clf;        
plot_mesh(Y,F, options);

exo3()

%% Insert your code here.
