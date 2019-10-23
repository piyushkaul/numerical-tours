
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/fastmarching_6_sampling_surf')

clear options;
name = 'bunny';
[vertex,faces] = read_mesh(name);
n = size(vertex,2);
options.name = name;

clf;
plot_mesh(vertex,faces, options);

landmarks = [100];

[D,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks);

clf; hold on;
options.face_vertex_color = mod( 20*D/max(D),1 );
plot_mesh(vertex,faces, options);
colormap jet(256);
h = plot3(vertex(1,landmarks), vertex(2,landmarks), vertex(3,landmarks), 'r.');
set(h, 'MarkerSize', 20);

[tmp,landmarks(end+1)] = max(D);

options.constraint_map = D;
[D1,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks,options);
D = min(D,D1);

clf; hold on;
options.face_vertex_color = mod( 20*D/max(D),1 );
plot_mesh(vertex,faces, options);
colormap jet(256);
h = plot3(vertex(1,landmarks), vertex(2,landmarks), vertex(3,landmarks), 'r.');
set(h, 'MarkerSize', 20);

exo1()

%% Insert your code here.

[D,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks);

[B,I,J] = unique(Q);
v = randperm(m)'; J = v(J);
clf; hold on;
options.face_vertex_color = J;
plot_mesh(vertex,faces, options);
colormap jet(256);
h = plot3(vertex(1,landmarks), vertex(2,landmarks), vertex(3,landmarks), 'k.');
set(h, 'MarkerSize', 15);

V = Q(faces); V = sort(V,1);
V = unique(V', 'rows')';
d = 1 + (V(1,:)~=V(2,:)) + (V(2,:)~=V(3,:));

I = find(d==3); I = sort(I);

z = zeros(n,1);
z(landmarks) = (1:m)';
facesV = z(V(:,I));

vertexV = vertex(:,landmarks);

options.method = 'slow';
options.verb = 0;
facesV = perform_faces_reorientation(vertexV,facesV, options);

clf;
options.face_vertex_color = [];
plot_mesh(vertexV,facesV, options);
shading faceted;

W = ones(n,1);
W(vertex(1,:)<median(vertex(1,:))) = .4;
options.W = W;

clf;
hold on;
options.face_vertex_color = W;
plot_mesh(vertex,faces, options);
colormap jet(256);

landmarks = [5000];
options.constraint_map = [];
[D,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks, options);

clf;
hold on;
options.face_vertex_color = mod( 20*D/max(D),1 );
plot_mesh(vertex,faces, options);
colormap jet(256);
h = plot3(vertex(1,landmarks), vertex(2,landmarks), vertex(3,landmarks), 'r.');
set(h, 'MarkerSize', 20);

exo2()

%% Insert your code here.

[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertex,faces,options);

C = abs(Cmin)+abs(Cmax);

clf;
hold on;
options.face_vertex_color = min(C,.1);
plot_mesh(vertex,faces, options);
colormap jet(256);

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.
