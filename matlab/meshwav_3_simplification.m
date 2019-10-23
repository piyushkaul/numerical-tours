
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/meshwav_3_simplification')

name = 'venus';
options.name = name;
[vertex,faces] = read_mesh(name);
n = size(vertex,2);

plot_mesh(vertex,faces,options);
shading faceted;

faces1 = faces;
vertex1 = vertex;

edges = compute_edges(faces1);
nedges = size(edges,2);

k = floor(rand*(nedges-1))+1;
e = edges(:,k);

vertex1(:,e(1)) = mean( vertex1(:,e),2 );
vertex1(:,e(2)) = Inf;

faces1(faces1==e(2)) = e(1);
a = sum( diff(sort(faces1))==0 );
faces1(:,a>0) = [];

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.
