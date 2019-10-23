
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/meshproc_4_fourier')

name = 'elephant-50kv';
[vertex,faces] = read_mesh(name);
options.name = name; 
n = size(vertex,2);

E = [faces([1 2],:) faces([2 3],:) faces([3 1],:)];
p = size(E,2);

W = sparse( E(1,:), E(2,:), ones(p,1) );
W = max(W,W');

D = spdiags(sum(W)', 0, n,n);
L = D-W;

nb = 80;
opts.disp = 0;
[U,S] = eigs(L,nb,'SM',opts);
S = diag(S);

[S,I] = sort(S, 'ascend');
U = real( U(:,I) );

clf;
plot(S); axis('tight');

ilist = round(linspace(3,nb, 6));
tau=2.2; % saturation for display
clf;
for i=1:length(ilist)
    v = real(U(:,ilist(i)));
    v = clamp( v/std(v),-tau,tau );
    options.face_vertex_color = v;
    subplot(2,3,i);
    plot_mesh(vertex,faces,options);
    shading interp; camlight; axis tight;
    colormap jet(256);
end

pvertex = vertex*U;

clf;
plot(pvertex'); axis('tight');
legend('X', 'Y', 'Z');

vertex1 = pvertex*U';

clf;
subplot(1,2,1);
plot_mesh(vertex,faces);
subplot(1,2,2);
plot_mesh(vertex1,faces);

exo1()

%% Insert your code here.

name = 'venus';
[vertex,faces] = read_mesh(name);
options.name = name;
n = size(vertex,2);

E = [faces([1 2],:) faces([2 3],:) faces([3 1],:)];
W = sparse( E(1,:), E(2,:), ones(size(E,2),1) );
W = max(W,W');
L = spdiags(sum(W)', 0, n,n) - W;

[U,S] = eig(full(L));
S = diag(S);
[S,I] = sort(S, 'ascend');
U = real( U(:,I) );

clf;
plot(S); axis('tight');

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

T = .05;

pvertexI = floor(abs(pvertex/T)).*sign(pvertex);

pvertexQ = sign(pvertexI) .* (abs(pvertexI)+.5) * T;

vertex1 = pvertexQ*U';

clf;
subplot(1,2,1);
plot_mesh(vertex,faces);
subplot(1,2,2);
plot_mesh(vertex1,faces);

t = min(pvertexI(:)):max(pvertexI(:));
h = hist( pvertexI(:), t );
h = max(h,1e-10); h = h/sum(h);

close; clf;
bar(t, h);
axis([-5 5 0 max(h)]);

E = -sum( log2(h).*h );

disp(['Nbr.bits per vertex = ' num2str(3*E,3)]);
disp(['Error,          SNR = ' num2str(snr(vertex,vertex1),3) 'dB']);

exo4()

%% Insert your code here.
