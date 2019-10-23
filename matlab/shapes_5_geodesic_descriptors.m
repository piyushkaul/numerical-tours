
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/shapes_5_geodesic_descriptors')

n = 200;
name = 'centaur1';
M = load_image(name,n);
M = perform_blurring(M,5);
M = double( rescale( M )>.5 );
if M(1)==1 
    M = 1-M;
end

clf;
imageplot(-M);

bound = compute_shape_boundary(M);
nbound = size(bound,2);

W = ones(n);
L = zeros(n)-Inf; L(M==1) = +Inf;

start_points = [95; 20];

options.constraint_map = [];
D0 = perform_fast_marching(W, start_points, options);
D0(M==0) = Inf;

clf;
options.display_levelsets = 1;
options.pstart = start_points;
options.nbr_levelsets = 30;
display_shape_function(D0, options);

options.constraint_map = L;
D = perform_fast_marching(W, start_points, options);

clf;
options.nbr_levelsets = 60;
display_shape_function(D, options);

exo1()

%% Insert your code here.

end_points = [27;112];
p = compute_geodesic(D,end_points);

ms = 30; lw = 3;
clf; hold on;
imageplot(1-M);
h = plot(end_points(2),end_points(1), '.b'); set(h, 'MarkerSize', ms);
h = plot(start_points(2),start_points(1), '.r'); set(h, 'MarkerSize', ms);
h = plot( p(2,:), p(1,:), 'g' ); set(h, 'LineWidth', lw);
axis ij;

exo2()

%% Insert your code here.

nb_samples = 600;
sel = round(linspace(1,nbound+1,nb_samples+1)); sel(end) = [];
samples = bound(:,sel);

exo3()

%% Insert your code here.

E = E/mean(E(:));

points = [[80;20] [95;112] [156;42]];
col = {'r', 'g', 'b', 'k'};
clf; hold on;
imageplot(-M);
for i=1:3
    h = plot(points(2,i), points(1,i), [col{i} '.']);
    set(h, 'MarkerSize', 40);
end
axis('ij');

clf;
col = {'r', 'g', 'b', 'k'};
for i=1:3
    subplot(3,1,i);    
    d = E(points(1,i),points(2,i), :); 
    u = hist(d(:), 15); axis tight;
    h = bar(u, col{i}); axis('tight');
    set(gca, 'XTickLabel', []);
end

clear A;
A{1} = max(E,[],3);
A{2} = min(E,[],3);
A{3} = mean(E,3);
A{4} = median(E,3);
titles = {'Max', 'Min', 'Mean', 'Median'};

nbr = [20 5 30 30];
options.pstart = [];
clf;
for i=1:4
    subplot(2,2,i);
    options.nbr_levelsets = nbr(i);
    display_shape_function(A{i}, options);
    title(titles{i});
end
colormap jet(256);

clf;
for i=1:4
    u = A{i}(M==1); u = u(u>0);
    subplot(4,1,i);
    hist(u, 40); axis('tight');
    title(titles{i});
end

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.
