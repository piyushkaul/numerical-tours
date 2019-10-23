
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/sparsity_4_dictionary_learning')

if not(exist('w'))
    w = 10;
end

n = w*w;

p = 2*n;

m = 20*p;

k = 4;

if not(exist('f'))
    f = rescale( crop(load_image('barb'),256) );
end
n0 = size(f,1);

clf;
imageplot(clamp(f));

q = 3*m;
x = floor( rand(1,1,q)*(n0-w) )+1;
y = floor( rand(1,1,q)*(n0-w) )+1;

[dY,dX] = meshgrid(0:w-1,0:w-1);
Xp = repmat(dX,[1 1 q]) + repmat(x, [w w 1]);
Yp = repmat(dY,[1 1 q]) + repmat(y, [w w 1]);
Y = f(Xp+(Yp-1)*n0);
Y = reshape(Y, [n q]);

Y = Y - repmat( mean(Y), [n 1] );

[tmp,I] = sort(sum(Y.^2), 'descend');
Y = Y(:,I(1:m));

ProjC = @(D)D ./ repmat( sqrt(sum(D.^2)), [w^2, 1] );
sel = randperm(m); sel = sel(1:p); 
D0 = ProjC( Y(:,sel) );
D = D0;

clf;
plot_dictionnary(D, [], [8 12]);

select = @(A,k)repmat(A(k,:), [size(A,1) 1]);
ProjX = @(X,k)X .* (abs(X) >= select(sort(abs(X), 'descend'),k));

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

clf;
plot_dictionnary(D,X, [8 12]);
