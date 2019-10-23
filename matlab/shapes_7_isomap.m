
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/shapes_7_isomap')

n = 1000;

x = rand(2,n);

v = 3*pi/2 * (.1 + 2*x(1,:));
X  = zeros(3,n);
X(2,:) = 20 * x(2,:);
X(1,:) = - cos( v ) .* v;
X(3,:) = sin( v ) .* v;

ms = 50;
lw = 1.5;
v1 = -15; v2 = 20;

clf;
scatter3(X(1,:),X(2,:),X(3,:),ms,v, 'filled'); 
colormap jet(256);
view(v1,v2); axis('equal'); axis('off');

D1 = repmat(sum(X.^2,1),n,1);
D1 = sqrt(D1 + D1' - 2*X'*X);

k = 6;

[DNN,NN] = sort(D1);
NN = NN(2:k+1,:);
DNN = DNN(2:k+1,:);

B = repmat(1:n, [k 1]);
A = sparse(B(:), NN(:), ones(k*n,1));

W = sparse(B(:), NN(:), DNN(:));

options.lw = lw;
options.ps = 0.01;
clf; hold on;
scatter3(X(1,:),X(2,:),X(3,:),ms,v, 'filled'); 
plot_graph(A, X, options);
colormap jet(256);
view(v1,v2); axis('equal'); axis('off');
zoom(.8);

D = full(W);
D = (D+D')/2;

D(D==0) = Inf;

D = D - diag(diag(D));

exo1()

%% Insert your code here.

Iremove = find(D(:,1)==Inf);

D(D==Inf) = 0;

exo2()

%% Insert your code here.

[U,L] = eig(Xstrain*Xstrain' / n);
Xstrain1 = U'*Xstrain;

Xstrain1(:,Iremove) = Inf;

clf; hold on;
scatter(Xstrain1(1,:),Xstrain1(2,:),ms,v, 'filled'); 
plot_graph(A, Xstrain1, options);
colormap jet(256);
axis('equal'); axis('off');

Y = cat(1, v, X(2,:));
Y(1,:) = rescale(Y(1,:), min(Xstrain(1,:)), max(Xstrain(1,:)));
Y(2,:) = rescale(Y(2,:), min(Xstrain(2,:)), max(Xstrain(2,:)));

clf; hold on;
scatter(Y(1,:),Y(2,:),ms,v, 'filled'); 
plot_graph(A,  Y, options);
colormap jet(256);
axis('equal'); axis('off'); 
camroll(90);

exo3()

%% Insert your code here.

clf;
plot(stress(1:end), '.-');
axis('tight');

[U,L] = eig(Xstress*Xstress' / n);
[L,I] = sort(diag(L));
U = U(:,I(2:3));

Xstress1 = U'*Xstress;

Xstress1(:,Iremove) = Inf;

clf; hold on;
scatter(Xstress1(1,:),Xstress1(2,:),ms,v, 'filled'); 
plot_graph(A, Xstress1, options);
colormap jet(256);
axis('equal'); axis('off');

exo4()

%% Insert your code here.
