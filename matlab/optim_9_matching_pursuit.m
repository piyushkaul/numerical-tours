
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/optim_9_matching_pursuit')

N = 200;
P = round(N/4);

Phi = randn(P,N);

Phi = Phi ./ repmat( sqrt(sum(Phi.^2)), [P 1] );

s = round(P/5);

I = randperm(N); I = I(1:s);
x0 = zeros(N,1); x0(I) = sign(randn(s,1));

sigma = 0.05 * norm(Phi*x0)/sqrt(P);

y = Phi*x0 + sigma*randn(P,1);

x = zeros(N,1);

c = Phi'*(y-Phi*x);

clf; 
stem( abs(c), 'b.' );
axis tight;

[~,i] = max(abs(c));

x(i) = x(i) + c(i);

exo1()

%% Insert your code here.

clf;
subplot(2,1,1);
x1 = x0; x1(x1==0) = NaN; stem(x1, 'b.'); axis tight;
subplot(2,1,2);
x1 = x; x1(x1==0) = NaN; stem(x1, 'b.'); axis tight;

clf;
I = find(x0~=0); J = setdiff(1:N,I);
clf; hold on;
h = plot(E, X(I,:)', '-'); set(h, 'LineWidth', 2);
h = plot(E, X(J,:)', 'k-'); set(h, 'LineWidth', 2);
axis tight; box on;

c = Phi'*(y-Phi*x);
[~,i] = max(abs(c));
x(i) = x(i) + c(i);

I = find(x~=0);
x(I) = pinv(Phi(:,I))*y;

exo2()

%% Insert your code here.

clf;
subplot(2,1,1);
x1 = x0; x1(x1==0) = NaN; stem(x1, 'b.'); axis tight;
subplot(2,1,2);
x1 = x; x1(x1==0) = NaN; stem(x1, 'b.'); axis tight;
