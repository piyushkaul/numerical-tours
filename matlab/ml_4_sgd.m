
addpath('toolbox_general')
addpath('solutions/ml_4_sgd')

f = {@(x)1/2*(x-1).^2, @(x)1/2*(x+1).^2};
F = @(x)f{1}(x)+f{2}(x);
df = {@(x)(x-1), @(x)(x+1)};

exo1()

%% Insert your code here.

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);
Xm = @(X)X-repmat(mean(X,1), [size(X,1) 1]);
Cov = @(X)Xm(X)'*Xm(X);

name = 'quantum';
load(['ml-' name]);
A = A(randperm(size(A,1)),:);
X = A(:,1:end-1);
y = A(:,end);

y = rescale(y,-1,1);

I = find(mean(abs(X))>1e-1); X = X(:,I);
X = X-repmat(mean(X),[size(X,1),1]);
X = X ./ repmat( sqrt(sum(X.^2)/size(X,1)), [size(X,1),1] );

[n,p] = size(X);

I = randperm(n); I = I(1:500);
options.disp_dim = 3;
clf; plot_multiclasses(X(I,:),y(I),options);

L = @(s,y)1/n * sum( log( 1 + exp(-s.*y) ) );
E = @(w,X,y)L(X*w,y);
theta = @(v)1 ./ (1+exp(-v));
nablaL = @(s,r)- 1/n * y.* theta(-s.*y);
nablaE = @(w,X,y)X'*nablaL(X*w,y);

exo2()

%% Insert your code here.

nablaEi = @(w,i)-y(i) .* X(i,:)' * theta( -y(i) * (X(i,:)*w)  );

l0 = 100;
tau0 = .05;

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.
