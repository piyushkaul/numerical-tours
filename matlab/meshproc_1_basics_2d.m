
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/meshproc_1_basics_2d')

n = 200;

vertex = 2*rand(2,n)-1;

faces = delaunay(vertex(1,:),vertex(2,:))';

clf;
subplot(1,2,1);
hh = plot(vertex(1,:),vertex(2,:), 'k.');
axis('equal'); axis('off');
set(hh,'MarkerSize',10);
title('Points');
subplot(1,2,2);
plot_mesh(vertex,faces);
title('Triangulation');

m = 20;
t = linspace(0,2*pi,m+1); t(end) = [];
vertexF = [cos(t);sin(t)];
vertex(:,1:m) = vertexF;
faces = delaunay(vertex(1,:),vertex(2,:))';

vertex1 = vertex;

faces1 = delaunay(vertex1(1,:),vertex1(2,:))';

E = [faces([1 2],:) faces([2 3],:) faces([3 1],:)];
p = size(E,2);

A = sparse( E(1,:), E(2,:), ones(p,1) );

d = 1./sum(A);
iD = spdiags(d(:), 0, n,n);
W = iD * A;

vertex1 = vertex1*W';

vertex1(:,1:m) = vertexF;

clf;
subplot(1,2,1);
plot_mesh(vertex,faces);
title('Before filering');
subplot(1,2,2);
plot_mesh(vertex1,faces1);
title('After filtering');

exo1()

%% Insert your code here.
