
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('toolbox_wavelet_meshes')
addpath('solutions/meshdeform_2_parameterization_sphere')

name = 'horse';
[vertex,faces] = read_mesh(name);
n = size(vertex,2);
m = size(faces,2);
clear options; options.name = name;

clf;
options.lighting = 1;
plot_mesh(vertex,faces,options);
shading faceted;

weight = 'conformal';
weight = 'combinatorial';
switch weight
    case 'conformal'
        W = make_sparse(n,n);
        for i=1:3
            i1 = mod(i-1,3)+1;
            i2 = mod(i  ,3)+1;
            i3 = mod(i+1,3)+1;
            pp = vertex(:,faces(i2,:)) - vertex(:,faces(i1,:));
            qq = vertex(:,faces(i3,:)) - vertex(:,faces(i1,:));
            % normalize the vectors
            pp = pp ./ repmat( sqrt(sum(pp.^2,1)), [3 1] );
            qq = qq ./ repmat( sqrt(sum(qq.^2,1)), [3 1] );
            % compute angles
            ang = acos(sum(pp.*qq,1));
            a = max(1 ./ tan(ang),1e-1); % this is *very* important
            W = W + make_sparse(faces(i2,:),faces(i3,:), a, n, n );
            W = W + make_sparse(faces(i3,:),faces(i2,:), a, n, n );
        end
    case 'combinatorial'
        E = [faces([1 2],:) faces([2 3],:) faces([3 1],:)];
        E = unique_rows([E E(2:-1:1,:)]')';
        W = make_sparse( E(1,:), E(2,:), ones(size(E,2),1) );
end

d = full( sum(W,1) );
D = spdiags(d(:), 0, n,n);
iD = spdiags(d(:).^(-1), 0, n,n);
tW = iD * W;

vertex1 = vertex;
vertex1 = vertex1 - repmat( mean(vertex1,2), [1 n] );
vertex1 = vertex1 ./ repmat( sqrt(sum(vertex1.^2,1)), [3 1] );

[normal,normalf] = compute_normal(vertex1,faces);

C = squeeze(mean(reshape(vertex1(:,faces),[3 3 m]), 2));

I = sum(C.*normalf);

disp(['Ratio of inverted triangles:' num2str(sum(I(:)<0)/m, 3) '%']);

options.name = 'none';
clf;
options.face_vertex_color = double(I(:)>0);
plot_mesh(vertex1,faces,options);
colormap gray(256); axis tight;
shading faceted;

vertex1 = vertex1*tW';
vertex1 = vertex1 ./ repmat( sqrt(sum(vertex1.^2,1)), [3 1] );

exo1()

%% Insert your code here.

clf;
plot(Edir/Edir(1));
axis('tight');

clf;
plot_mesh(vertex1,faces);
colormap gray(256); axis tight;
shading faceted;

clf;
plot(ninvert/m); axis tight;

q = 128;
options.verb = 0;
M = perform_sgim_sampling(vertex, vertex1, faces, q, options);

clf;
plot_geometry_image(M, 1, 1);
axis equal;
colormap gray(256);
view(134,-61);

exo2()

%% Insert your code here.

vertex1 = vertex; 
vertex1 = vertex1 - repmat( mean(vertex1,2), [1 n] );
vertex1 = vertex1 ./ repmat( sqrt(sum(vertex1.^2,1)), [3 1] );

eta = .5;

A = squeeze(mean(reshape(vertex1(:,faces),[3 3 m]), 2));

E = zeros(1,m);
for i=1:3
    i1 = mod(i,3)+1;
    % directed edge
    u = vertex(:,faces(i,:)) - vertex(:,faces(i1,:));
    % norm squared
    u = sum(u.^2);
    % weights between the vertices
    w = W(faces(i,:) + (faces(i1,:)-1)*n);
    E = E + w.*u;
end

G = zeros(3,n);
for j=1:m
    f = faces(:,j);
    Alpha = A(:,j);
    alpha = norm(Alpha);
    for i=1:3
        i1 = mod(i  ,3)+1;
        i2 = mod(i+1,3)+1;
        % directed edges
        u1 = vertex(:,f(i)) - vertex(:,f(i1));
        u2 = vertex(:,f(i)) - vertex(:,f(i2));
        % weights between the vertices
        w1 = W(f(i) + (f(i1)-1)*n);
        w2 = W(f(i) + (f(i2)-1)*n);
        G(:,f(i)) = G(:,f(i)) + (w1*u1 + w2*u2) ./ alpha^2 - Alpha/alpha^4 * E(j);
    end
end

vertex1 = vertex1 - eta*G;
vertex1 = vertex1 ./ repmat( sqrt(sum(vertex1.^2,1)), [3 1] );

exo3()

%% Insert your code here.

clf;
plot(Edir/Edir(1));
