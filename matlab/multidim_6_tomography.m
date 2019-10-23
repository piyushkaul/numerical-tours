
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/multidim_6_tomography')

name = 'vessels';
options.nbdims = 3;
M = read_bin(name, options);
M = rescale(M);

n = size(M,1);

M = M(1:2:n,1:2:n,1:2:n);
n = n/2;

slices = round(linspace(10,n-10,4));
clf;
for i=1:length(slices)
    s = slices(i);
    imageplot( M(:,:,s), strcat(['Z=' num2str(s)]), 2,2,i );
end

sel = 1:2:n;
clf;
isosurface( M(sel,sel,sel), .5);
axis('off');

P = 12;

if P==12
    tau = 0.8506508084;
    one = 0.5257311121;
    S = [  tau,  one,    0;
        -tau,  one,    0
        -tau, -one,    0;
        tau, -one,    0;
        one,   0 ,  tau;
        one,   0 , -tau;
        -one,   0 , -tau;
        -one,   0 ,  tau;
        0 ,  tau,  one;
        0 , -tau,  one;
        0 , -tau, -one;
        0 ,  tau, -one ]';
else
    S = randn(3,P);
    S = S ./ repmat( sqrt( sum(S.^2,1) ), [3 1]);
end

clf;
plot3(S(1,:), S(2,:), S(3,:), '.');
axis equal;

x = [0:n/2-1, -n/2:-1];
[X,Y,Z] = ndgrid(x,x,x);
mask = zeros(n,n,n);
epsilon = .5;
for i=1:P
    q = S(:,i);
    d = q(1)*X + q(2)*Y + q(3)*Z;
    mask( abs(d)<=epsilon ) = 1;
end

F = fftn(M);
y = F(mask==1);

Q = length(y);
disp(strcat(['Number of measurements Q=' num2str(Q) '.']));
disp(strcat(['Sub-sampling Q/N=' num2str(length(y)/n^3,2) '.']));

F1 = zeros(n,n,n);
F1(mask==1) = y;
M1 = real( ifftn(F1) );

clf;
for i=1:length(slices)
    s = slices(i);
    imageplot( clamp(M1(:,:,s)), strcat(['Z=' num2str(s)]), 2,2,i );
end

disp(['Pseudo-inverse, SNR=' num2str(snr(M,M1),4) 'dB.']);

epsilon = 1e-2;

G = grad(M);

clf;
for i=1:length(slices)
    s = slices(i);
    imageplot( squeeze(G(:,:,s,:)), strcat(['Z=' num2str(s)]), 2,2,i );
end

d = sqrt(sum(G.^2,4)+epsilon^2);

clf;
for i=1:length(slices)
    s = slices(i);
    imageplot( squeeze(d(:,:,s,:)), strcat(['Z=' num2str(s)]), 2,2,i );
end

tv = sum(d(:));
disp(['TV norm=' num2str(tv) '.']);

tau = epsilon*.2;

Mtv = M1;

G = grad(Mtv);
d = sqrt(sum(G.^2,4)+epsilon^2);
dG = -div( G ./ repmat(d, [1 1 1 3]) );

clf;
for i=1:length(slices)
    s = slices(i);
    imageplot( dG(:,:,s), strcat(['Z=' num2str(s)]), 2,2,i );
end

Mtv = Mtv - tau*dG;

F = fftn(Mtv);
F(mask==1) = y;
Mtv = real(ifftn(F));

exo1()

%% Insert your code here.

clf;
for i=1:length(slices)
    s = slices(i);
    imageplot( clamp(Mtv(:,:,s)), strcat(['Z=' num2str(s)]), 2,2,i );
end

disp(['Total variation, SNR=' num2str(snr(M,Mtv),4) 'dB.']);
