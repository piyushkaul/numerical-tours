
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/fastmarching_9_heuristics')

n = 300;
name = 'road2';
M = rescale(load_image(name,n));

clf;
imageplot(M);

pstart = [14;161];
pend = [293;148];

W = abs(M-M(pstart(1),pstart(2)));
W = rescale(W, 1e-2,1);

clf;
imageplot(W);

[D,S] = perform_fast_marching(1./W, pstart);

p = compute_geodesic(D,pend);

clf; hold on;
imageplot(convert_distance_color(D,M), 'Distance');
h = plot(p(2,:),p(1,:), '.k'); set(h, 'LineWidth', 2);
h = plot(pstart(2),pstart(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(pend(2),pend(1), '.b'); set(h, 'MarkerSize', 25);
axis ij;

[H,S] = perform_fast_marching(1./W, pend);

clf; hold on;
imageplot(convert_distance_color(H,M), 'Distance');
h = plot(p(2,:),p(1,:), '.k'); set(h, 'LineWidth', 2);
h = plot(pstart(2),pstart(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(pend(2),pend(1), '.b'); set(h, 'MarkerSize', 25);
axis ij;

exo1()

%% Insert your code here.

weight = .9;
options.end_points = pend;
options.heuristic = weight*H;
options.nb_iter_max = Inf;
options.constraint_map = Inf+zeros(n);
[D,S] = perform_fast_marching(1./W, pstart, options);

I = find(S<0);
U = cat(3,M,M,M);
U(I) = 1; U([I+n^2, I+2*n^2]) = U([I+n^2, I+2*n^2])*.3;
clf; hold on;
imageplot(U);
h = plot(p(2,:),p(1,:), '.k'); set(h, 'LineWidth', 2);
h = plot(pstart(2),pstart(1), '.g'); set(h, 'MarkerSize', 25);
h = plot(pend(2),pend(1), '.b'); set(h, 'MarkerSize', 25);
axis ij;

exo2()

%% Insert your code here.

q = 10;
landmarks = floor(rand(2,q)*n)+1;

Dland = zeros(n,n,q);
for i=1:q
    Dland(:,:,i) = perform_fast_marching(1./W, landmarks(:,i));
end

Dend = Dland( pend(1), pend(2), :);
H = max(abs(Dland-repmat(Dend, [n n 1])), [], 3);

clf;
hold on;
imageplot(H);
contour(H, 10, 'k', 'LineWidth', 2);
colormap jet(256);
h = plot(landmarks(1,:), landmarks(2,:), 'y.');
set(h, 'MarkerSize', 15);
axis ij;

[H0,S] = perform_fast_marching(1./W, pend);
clf;
hold on;
imageplot(H0);
contour(H0, 10, 'k', 'LineWidth', 2);
colormap jet(256);
axis ij;

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.

exo6()

%% Insert your code here.
