
addpath('toolbox_general')
addpath('solutions/ml_5_nn')

dotp = @(x,y)sum(x(:).*y(:));
max2 = @(S)repmat(max(S,[],2), [1 size(S,2)]);
SM = @(S)exp(S) ./ repmat( sum(exp(S),2), [1 size(S,2)]);
SM = @(S)SM(S-max2(S));

n0 = 100; % number of points per class
p = 2;   % dimensionality
k = 3;   % number of classes
n = n0*k; % Total number of points
x = zeros(p,n);
y = zeros(1,n);
for j=1:k
    I = n0*(j-1)+1:n0*j;
    r = linspace(0.0,1,n0); % radius
    t = linspace(j*4,(j+1)*4,n0) + randn(1,n0)*0.2; % angle
    x(:,I) = [r.*sin(t); r.*cos(t)];
    y(1,I) = j;
end

col = {'r' 'g' 'b'};
clf;  hold on;
for j=1:k
    I = find(y==j);
    plot(x(1,I), x(2,I), '.', 'color', col{j}, 'MarkerSize', 20);
end
axis equal; axis off;

Y = double( repmat((1:k)', [1,n]) == repmat(y, [k,1]) );

dotp = @(x,y)sum(x(:).*y(:));
max2 = @(S)repmat(max(S,[],2), [1 size(S,2)]);

LSE = @(S)log( sum(exp(S), 2) );
LSE = @(S)LSE( S-max2(S) ) + max(S,[],2);

SM = @(S)exp(S) ./ repmat( sum(exp(S),2), [1 size(S,2)]);
SM = @(S)SM(S-max2(S));

loss.F = @(z,y)sum(LSE(z')) - dotp(z,y);

loss.G = @(z,y)SM(z')' -  y;

rho.F  = @(u)atan(u);
rho.G = @(u)1./(1+u.^2);

t = linspace(-5,5,201);
clf;
plot(t, rho.F(t), 'LineWidth', 2);
axis tight;

D = [p 15 k]; % here a single hidden layer

R = length(D)-1; 
A = {}; b = {}; 
for r=1:R
    A{r} = randn(D(r+1),D(r));
    b{r} = randn(D(r+1),1);
end

X = {};
X{1} = x;
for r=1:R
    X{r+1} = rho.F( A{r}*X{r}+repmat(b{r},[1 n]) );
end
L = loss.F(X{R+1},Y);

gx = loss.G(X{R+1},Y);

gA = {}; gb = {};
for r=R:-1:1
    M = rho.G(A{r}*X{r}+repmat(b{r},[1 n])) .* gx;
    % nabla_X{r} 
    gx = A{r}' * M;
    % nabla_A{r}
    gA{r} =  M * X{r}';
    % nabla_b{r}
    gb{r} =  sum(M,2);
end

exo1()

%% Insert your code here.

q = 100;
t = linspace(-1,1,q);
[Yg,Xg] = meshgrid(t,t);
Z = [Xg(:)';Yg(:)'];

V = EvalNN(Z,[], A,b,loss,rho); 
U = reshape(SM(V{end}'), [q q k]);

col = [ [1 0 0]; [0 1 0]; [0 0 1]; [0 0 0]; [0 1 1]; [1 0 1]; [1 1 0]; ...
    [1 .5 .5]; [.5 1 .5]; [.5 .5 1]  ]';
R = zeros(q,q,3);
for i=1:k
    for a=1:3
        R(:,:,a) = R(:,:,a) + U(:,:,i) .* col(a,i);
    end
end

clf;  hold on;
imagesc(t,t,permute(R, [2 1 3]));
for j=1:k
    I = find(y==j);
    plot(x(1,I), x(2,I), '.', 'color', col(:,j)*.8, 'MarkerSize', 20);
end
axis equal; axis off;

exo2()

%% Insert your code here.
