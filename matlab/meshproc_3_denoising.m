
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/meshproc_3_denoising')

clear options;
name = 'nefertiti';
name = 'elephant-50kv';
options.name = name; % useful for displaying
[X0,F] = read_mesh(name);

n = size(X0,2);
m = size(F,2);

options.lighting = 1;
% clf;
plot_mesh(X0,F, options);

rho = 0.015;

N = compute_normal(X0,F);

X = X0 + repmat(rho*randn(1,n),[3,1]).*N;

clf;
plot_mesh(X,F,options); axis('tight');

E = [F([1 2],:) F([2 3],:) F([3 1],:)];

E = unique_rows([E E(2:-1:1,:)]')';

E0 = E(:,E(1,:)<E(2,:));

p0 = size(E0,2);

[tmp,I] = sort(E0(1,:)');
E0 = E0(:,I);

disp(['#vertices=' num2str(n) ', #faces=' num2str(m) ', #edges=' num2str(p0) '.']);

W = make_sparse( E(1,:), E(2,:), ones(size(E,2),1) );

d = full( sum(W,1) );

clf;
hist(d,min(d):max(d)); 
axis('tight');

D = spdiags(d(:), 0, n,n);
iD = spdiags(d(:).^(-1), 0, n,n);

tW = iD * W;

L = D - W;

G = make_sparse( [1:p0 1:p0], [E0(1,:) E0(2,:)], [ones(1,p0) -ones(1,p0)] );

clf;
subplot(1,2,1);
spy(W); title('W');
subplot(1,2,2);
spy(G); title('G');

err = norm( G'*G - L, 'fro');
disp(['Factorization error (should be 0) = ' num2str(err,2) '.']);

M = load_image('lena',256);

v = X0 - repmat(mean(X0,2), [1 n]);
theta = acos(v(1,:)./sqrt(sum(v.^2)))/pi;
phi = (atan2(v(2,:),v(3,:))/pi+1)/2;

x = linspace(0,1,size(M,1));
f = interp2(x,x,M',theta,phi)';

options.face_vertex_color = f(:);
clf;
plot_mesh(X0,F, options);
lighting none;

exo1()

%% Insert your code here.

niter = 5;
X1 = X;
for i=1:niter
    X1 = X1*tW'; 
end

pnoisy = snr(X0,X);
pfilt  = snr(X0,X1);
disp(strcat(['Noisy=' num2str(pnoisy,2) 'dB, denoised=' num2str(pfilt,2) 'dB.']));

clf;
plot_mesh(X1,F, options);
axis('tight'); shading('interp');

exo2()

%% Insert your code here.

clf;
plot(0:length(err)-1, err, '.-'); axis('tight');
set_label('Iteration', 'SNR');

tL = iD * L;

tau = .2;

Tmax = 40;

niter = ceil(Tmax/tau);

Xt = X;

Xt = Xt - tau*Xt*tL';

exo3()

%% Insert your code here.

t = linspace(0,Tmax,niter);
clf;
plot(t, err); axis('tight');
set_label('Time', 'SNR');

mu = 10;

A = speye(n,n)+mu*L;

Xmu = X;
for i=1:3
    b = X(i,:)';
    Xmu(i,:) = perform_cg(A,b)';
end

clf;
plot_mesh(Xmu,F, options);

exo4()

%% Insert your code here.
