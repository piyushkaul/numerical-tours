
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('toolbox_additional')
addpath('solutions/meshproc_6_volumetric')

[vertex,faces] = read_tet('hand.tet');

clear options;
options.plot_points = 1;
clf; plot_mesh(vertex,faces,options);

options.cutting_plane = [0 0 1];
options.plot_points = 0;
clf; plot_mesh(vertex,faces,options);

options.cutting_plane = [0 -1 0];
options.plot_points = 0;
options.cutting_offs = -.2;
options.face_vertex_color = vertex(1,:)';
clf; plot_mesh(vertex,faces,options);
view(-20,45);zoom(.8);
colormap jet(256);

exo1()

%% Insert your code here.

n = size(W,1);
D = spdiags(sum(W)', 0, n,n);
L = D-W;

i = round( rand(20,1)*n ) + 1;
b = zeros(n,1); 
b(i) = (-1).^(1:length(i));

L1 = L;
L1(i,:) = 0; L1(i+(i-1)*n) = 1;
v = L1\b;

options.face_vertex_color = v;
clf; plot_mesh(vertex,faces,options);
view(-20,45);zoom(.8);
shading interp;
colormap jet(256);
