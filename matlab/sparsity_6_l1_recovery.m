
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/sparsity_6_l1_recovery')

remove_diag = @(C)C-diag(diag(C));
Correlation = @(Phi)remove_diag(abs(Phi'*Phi));

maxall = @(C)max(C(:));
mu = @(Phi)maxall(Correlation(Phi));

Coh = @(Phi,k)(k * mu(Phi)) / ( 1 - (k-1) * mu(Phi) );

normalize = @(Phi) Phi ./ repmat(sqrt(sum(Phi.^2)), [size(Phi,1) 1]);
PhiRand = @(P,N)normalize(randn(P,N));
Phi = PhiRand(250,1000);

c = mu(Phi);
fprintf('Coherence:    %.2f\n', c);
fprintf('Sparsity max: %d\n', floor(1/2*(1+1/c)) );

exo1()

%% Insert your code here.

supp   = @(s)find(abs(s)>1e-5);
cosupp = @(s)find(abs(s)<1e-5);

PsiI = @(Phi,I)Phi(:, setdiff(1:size(Phi,2),I) )' * pinv(Phi(:,I))';

F = @(Phi,s)norm(PsiI(Phi,supp(s))*s(supp(s)), 'inf');

erc = @(Phi,I)norm(PsiI(Phi,I), 'inf');

g = @(C,I)sum(C(:,I),2);
werc_g = @(g,I,J)max(g(J)) / (1-max(g(I)));
werc = @(Phi,I)werc_g( g(Correlation(Phi),I), I, setdiff(1:size(Phi,2),I) );

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

minmax = @(v)deal(1-min(v),max(v)-1);
ric = @(A)minmax(eig(A'*A));

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.

sigma = 6;
g = @(x)(1-x.^2/sigma^2).*exp(-x.^2/(2*sigma^2));

P = 1024;
[Y,X] = meshgrid(1:P,1:P);
Phi = normalize(g(mod(X-Y+P/2, P)-P/2));

eta = 2;
N = P/eta;
Phi = Phi(:,1:eta:end);

c = Phi'*Phi;
c = abs(c(:,end/2));
clf;
h = plot(c(end/2-50:end/2+50), '.-'); set(h, 'LineWidth', 1);
axis tight;

twosparse = @(d)circshift([1; zeros(d,1); -1; zeros(N-d-2,1)], round(N/2-d/2));

x0 = twosparse(50);
clf; 
subplot(2,1,1);
h = plot(x0, 'r'); axis tight;
subplot(2,1,2);
h = plot(Phi*x0, 'b'); axis tight;
set(h, 'LineWidth', 2);

exo6()

%% Insert your code here.
