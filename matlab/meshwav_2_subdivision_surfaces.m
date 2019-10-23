
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('toolbox_wavelet_meshes')
addpath('solutions/meshwav_2_subdivision_surfaces')

[vertex1,face1] = compute_base_mesh('oct');
[vertex0,face0] = compute_base_mesh('ico');

clf;
subplot(1,2,1);
plot_mesh(vertex1,face1); 
shading('faceted'); lighting('flat'); view(3); axis('tight');
subplot(1,2,2);
plot_mesh(vertex0,face0); 
shading('faceted'); lighting('flat'); view(3); axis('tight');

face = face0;
vertex = vertex0;

edge = compute_edges(face);

n = size(vertex,2); 
ne = size(edge,2);

A = sparse([edge(1,:);edge(2,:)],[edge(2,:);edge(1,:)],[n+(1:ne);n+(1:ne)],n,n);
v12 = full( A( face(1,:) + (face(2,:)-1)*n ) );
v23 = full( A( face(2,:) + (face(3,:)-1)*n ) );
v31 = full( A( face(3,:) + (face(1,:)-1)*n ) );

face = [  cat(1,face(1,:),v12,v31),...
    cat(1,face(2,:),v23,v12),...
    cat(1,face(3,:),v31,v23),...
    cat(1,v12,v23,v31)   ];

vertex = [vertex, (vertex(:,edge(1,:))+vertex(:,edge(2,:)))/2 ];

d = sqrt( sum(vertex.^2,1) );
vertex = vertex ./ repmat( d, [size(vertex,1) 1]);

clf;
subplot(1,2,1);
plot_mesh(vertex0,face0); 
shading('faceted'); lighting('flat'); view(3); axis('tight');
subplot(1,2,2);
plot_mesh(vertex,face); 
shading('faceted'); lighting('flat'); view(3); axis('tight');

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

name = 'mannequin';
[vertex0,face0] = read_mesh(name);

options.name = name;
clf;
plot_mesh(vertex0,face0,options);
shading('faceted'); lighting('flat'); axis('tight');

face = face0;
vertex = vertex0;

edge = compute_edges(face);
n = size(vertex,2);
ne = size(edge,2);

A = sparse([edge(1,:);edge(2,:)],[edge(2,:);edge(1,:)],[n+(1:ne);n+(1:ne)],n,n);
v12 = full( A( face(1,:) + (face(2,:)-1)*n ) );
v23 = full( A( face(2,:) + (face(3,:)-1)*n ) );
v31 = full( A( face(3,:) + (face(1,:)-1)*n ) );

face_old = face;
face = [  cat(1,face(1,:),v12,v31),...
    cat(1,face(2,:),v23,v12),...
    cat(1,face(3,:),v31,v23),...
    cat(1,v12,v23,v31)   ];

global vring e2f fring facej;
vring = compute_vertex_ring(face);
e2f = compute_edge_face_ring(face_old);
fring = compute_face_ring(face_old);
facej = face_old;

for k=n+1:n+ne
    [e,v,g] = compute_butterfly_neighbors(k, n);
    vertex(:,k) = 1/2*sum(vertex(:,e),2) + 1/8*sum(vertex(:,v),2) - 1/16*sum(vertex(:,g),2);
end

clf;
subplot(1,2,1);
plot_mesh(vertex0,face0,options); 
shading('faceted'); lighting('flat'); axis('tight');
subplot(1,2,2);
plot_mesh(vertex,face,options); 
shading('faceted'); lighting('flat'); axis('tight');

exo3()

%% Insert your code here.

clf;
plot_mesh(vertex,face,options); 
shading('interp'); lighting('phong'); axis('tight');

clf;
plot_mesh(vertex,face,options); 
shading('faceted'); lighting('phong'); axis('tight');

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.

clf;
plot_mesh(vertex,face,options); 
shading('interp'); lighting('phong'); axis('tight');

clf;
plot_mesh(vertex,face,options); 
shading('faceted'); lighting('phong'); axis('tight');

exo6()

%% Insert your code here.
