
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/optimaltransp_4_matching_sliced')

options.rows=1;
P = @(f,g)perform_hist_eq(f,g,options);

N = 300;

d = 2;

f = rand(2,N)-.5;

theta = 2*pi*rand(1,N);
r = .8 + .2*rand(1,N);
g = [cos(theta).*r; sin(theta).*r];

plotp = @(x,col)plot(x(1,:)', x(2,:)', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', 2);

clf; hold on;
plotp(f, 'b');
plotp(g, 'r');
axis('off'); axis('equal');

tau = .2;

f1 = f;

[Theta,~] = qr(randn(d));

f1 = (1-tau)*f1 + tau * Theta * P(Theta'*f1, Theta'*g);

clf; hold on;
plotp(f, 'b');
plotp(g, 'r');
plotp(f1, 'y');
axis('off'); axis('equal');

exo1()

%% Insert your code here.

clf;
hold on;
h = plot([f(1,:);f1(1,:)], [f(2,:);f1(2,:)], 'k');
set(h, 'LineWidth', 2);
plotp(f, 'b');
plotp(g, 'r');
axis('off'); axis('equal');

t = .5;
ft = (1-t)*f + t*f1;

clf;
hold on;
plotp(f, 'b');
plotp(g, 'r');
plotp(ft, [t 0 1-t]);
axis('off'); axis('equal');

exo2()

%% Insert your code here.

n = 128;
N = n*n;
d = 3;

F = rescale( load_image('hibiscus', n) );
G = rescale( load_image('flowers', n) );

clf;
imageplot({F G});

f = reshape(F, [n*n 3])';
g = reshape(G, [n*n 3])';

quantize = @(A,Q)1+round((Q-1)*A);
J = @(I,Q)I(1,:)' + Q*(I(2,:)'-1);
hist2d = @(f,Q)reshape( accumarray(J(quantize(f,Q),Q), ones(1,N), [Q*Q 1], @sum), [Q Q]);

Q = 60;

func = @(a)log(a+3);
clf;
imageplot({ func(hist2d(f(1:2,:),Q)), func(hist2d(g(1:2,:),Q)) });

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.

clf;
imageplot(F1);

exo5()

%% Insert your code here.
