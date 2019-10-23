
addpath('toolbox_general')
addpath('solutions/ml_2_regression')

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);
Xm = @(X)X-repmat(mean(X,1), [size(X,1) 1]);
Cov = @(X)Xm(X)'*Xm(X);

name = 'prostate';
load(['ml-' name]);

A = A(randperm(size(A,1)),:);

X = A(:,1:end-2);
y = A(:,end-1);
c = A(:,end);

[n,p] = size(X);
fprintf('n=%d, p=%d\n', n,p);

I0 = find(c==1); % train
I1 = find(c==0); % test
n0 = length(I0); n1 = n-n0;
X0 = X(I0,:); y0 = y(I0);
X1 = X(I1,:); y1 = y(I1);

RepSub = @(X,u)X-repmat(u, [size(X,1) 1]);
RepDiv = @(X,u)X./repmat(u, [size(X,1) 1]);
mX0 = mean(X0); sX0 = std(X0);
X0 = RepSub(X0,mX0);  X1 = RepSub(X1,mX0);
X0 = RepDiv(X0,sX0);  X1 = RepDiv(X1,sX0);

m0 = mean(y0);
y0 = y0-m0;
y1 = y1-m0;

clf;
imagesc(Cov(X0));

clf;
bar(X0'*y0);
axis tight;
SetAR(1/2);

[U,D,V] = svd(Xm(X0),'econ');
Z = Xm(X0) * V;

clf; 
plot(diag(D), '.-', 'LineWidth', 2, 'MarkerSize', 30);
axis tight;
SetAR(1/2);

pmax = min(p,8);
k = 0;
clf;
for i=1:pmax
    for j=1:pmax
        k = k+1;
        subplot(pmax,pmax,k);
        if i==j
            hist(X0(:,i),6);
            axis tight;
        else
            plot(X0(:,j),X0(:,i), '.');
            axis tight;
        end
        set(gca, 'XTick', [], 'YTick', [] );
        axis tight;
        if i==1
            title(class_names{j});
        end
    end
end

options.disp_dim = 3;
clf; plot_multiclasses(X,ones(n,1),options);
SetAR(1);

col = {'b' 'g' 'r' 'c' 'm' 'y' 'k'};
clf;
for i=1:min(p,3)
    subplot(3,1,i);
    plot(Z(:,i), y0, '.', 'Color', col{i}, 'MarkerSize', 20);
    axis tight;
end

w = (X0'*X0) \ (X0'*y0);

clf;
plot( [y1 X1*w], '.-', 'MarkerSize', 20);
axis tight;
legend('y', 'X_1 w');

E = norm(X1*w-y1) / norm(y1);
fprintf('Relative prediction error: %.3f\n', E);

lambda = .2*norm(X0)^2;
w = (X0'*X0+lambda*eye(p)) \ (X0'*y0);
w1 = X0'*( (X0*X0'+lambda*eye(n0)) \ y0 );
fprintf('Error (should be 0): %.4f\n', norm(w-w1)/norm(w));

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

J = @(w,lambda)1/2*norm(X0*w-y0)^2 + lambda*norm(w,1);

Soft = @(x,s)max(abs(x)-s,0).*sign(x);

t = linspace(-5,5,201);
clf; plot(t,Soft(t,2), 'LineWidth', 2); 
axis tight;
SetAR(1/2);

tau = 1.5/norm(X0)^2;

lambda = max(abs(X0'*y0))/10;

w = zeros(p,1);

C = X0'*X0;
u = X0'*y0;
ISTA = @(w,lambda,tau)Soft( w-tau*( C*w-u ), lambda*tau );
w = ISTA(w,lambda,tau);

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.

exo6()

%% Insert your code here.

B = 3;
n = 500; p = 2;
X = 2*B*rand(n,2)-B;
rho = .5; % noise level
y = peaks(X(:,1), X(:,2)) + randn(n,1)*rho;

clf;
scatter(X(:,1), X(:,2), ones(n,1)*20, y, 'filled');
colormap jet(256); 
axis equal; axis([-B B -B B]); box on;

distmat = @(X,Z)bsxfun(@plus,dot(X',X',1)',dot(Z',Z',1))-2*(X*Z');

sigma = .3;
kappa = @(X,Z)exp( -distmat(X,Z)/(2*sigma^2) );

K = kappa(X,X);

lambda = 0.01;
h = (K+lambda*eye(n))\y;

Y = @(x)kappa(x,X)*h;

q = 101;
t = linspace(-B,B,q);
[v,u] = meshgrid(t,t);
Xn = [u(:), v(:)];

yn = reshape(Y(Xn),[q,q]);
clf;
imagesc(t,t,yn); axis image; axis off; 
colormap jet(256);

exo7()

%% Insert your code here.

exo8()

%% Insert your code here.
