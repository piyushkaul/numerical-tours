
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/optim_8_homotopy')

N = 200;
P = round(N/4);

Phi = randn(P,N);

s = round(P/5);

I = randperm(N); I = I(1:s);
x0 = zeros(N,1); x0(I) = sign(randn(s,1));

sigma = 0.05 * norm(Phi*x0)/sqrt(P);

y = Phi*x0 + sigma*randn(P,1);

C = Phi'*y;
[lambda,I] = max(abs(C));
x = zeros(N,1);

J = setdiff(1:N, I);

c = Phi'*(y-Phi*x);

clf; hold on;
stem( I, c(I)/lambda, 'b.' );
stem( J, c(J)/lambda, 'r.' );
plot([1 N], [1 1], 'k--');
plot([1 N],-[1 1], 'k--');
axis([1 N -1.05 1.05]);

sI = sign(c(I));

d = zeros(N,1);
d(I) = (Phi(:,I)'*Phi(:,I)) \ sI;

v = Phi(:,I)*d(I);

w = ( lambda-c(J) ) ./ ( 1 - Phi(:,J)'*v );
gamma1 = min(w(w>0));
if not(isempty(gamma1))
    i1 = J( w==gamma1 );
end

w = ( lambda+c(J) ) ./ ( 1 + Phi(:,J)'*v );
gamma2 = min(w(w>0));
if not(isempty(gamma2))
    i2 = J( w==gamma2 );
end

w = -x(I)./d(I);
gamma3 = min(w(w>0));
if not(isempty(gamma3))
    i3 = I( w==gamma3 );
end

gamma = min([gamma1 gamma2 gamma3]);

x = x + gamma*d;

lambda = lambda - gamma;

if gamma==gamma1
    I = [I i1];
elseif gamma==gamma2
    I = [I i2];
elseif gamma==gamma3
    I(I==i3) = [];
end

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
h = plot(Lambda, X(I,:)', '-'); set(h, 'LineWidth', 2);
h = plot(Lambda, X(J,:)', 'k-'); set(h, 'LineWidth', 2);
axis tight;

clf;
plot(Lambda, Sparsity, '.-');
axis tight;

exo2()

%% Insert your code here.
