
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/meshdeform_5_deformation')

p = 150;
[Y,X] = meshgrid(linspace(-1,1,p),linspace(-1,1,p));
vertex0 = [X(:)'; Y(:)'; zeros(1,p^2)];
n = size(vertex0,2);

I = reshape(1:p^2,p,p);
a = I(1:p-1,1:p-1); b = I(2:p,1:p-1); c = I(1:p-1,2:p);
d = I(2:p,1:p-1); e = I(2:p,2:p); f = I(1:p-1,2:p);
faces = cat(1, [a(:) b(:) c(:)], [d(:) e(:) f(:)])';

sigma = .03;
h = .35;
q = 8;

t = linspace(-1,1,q+2); t([1 length(t)]) = [];
vertex = vertex0;
for i=1:q
    for j=1:q
        d = (X(:)'-t(i)).^2 + (Y(:)'-t(j)).^2;
        vertex(3,:) = vertex(3,:) + h * exp( -d/(2*sigma^2)  );
    end
end

clf;
plot_mesh(vertex,faces);
view(3);

W = sparse(n,n);
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
   u = cot(ang);
   u = clamp(u, 0.01,100);
   W = W + sparse(faces(i2,:),faces(i3,:),u,n,n);
   W = W + sparse(faces(i3,:),faces(i2,:),u,n,n);
end

d = full( sum(W,1) );
D = spdiags(d(:), 0, n,n);
L = D - W;

I = find( abs(X(:))==1 | abs(Y(:))==1 );

Delta0 = zeros(3,n);
d = ( vertex(1,I) + vertex(2,I) ) / 2;
Delta0(3,I) = sign(d) .* abs(d).^3;

L1 = L;
L1(I,:) = 0;
L1(I + (I-1)*n) = 1;

Delta = ( L1 \ Delta0' )';

vertex1 = vertex+Delta;

clf;
plot_mesh(vertex1,faces);
view(-100,15);

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.

normal = compute_normal(vertex0,faces);
d = repmat( sum(normal .* (vertex-vertex0)), [3 1]);

exo5()

%% Insert your code here.

exo6()

%% Insert your code here.

exo7()

%% Insert your code here.

exo8()

%% Insert your code here.

exo9()

%% Insert your code here.
