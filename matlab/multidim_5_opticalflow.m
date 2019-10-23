
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/multidim_5_opticalflow')

n = 256;
name = 'lena';
M1 = rescale( load_image(name,n) );

theta = .03 * pi/2;

[Y,X] = meshgrid(1:n,1:n);

X1 = (X-n/2)*cos(theta) + (Y-n/2)*sin(theta) + n/2;
Y1 =-(X-n/2)*sin(theta) + (Y-n/2)*cos(theta) + n/2;

X1 = mod(X1-1,n)+1;
Y1 = mod(Y1-1,n)+1;

M2 = interp2(Y,X,M1,Y1,X1);
M2(isnan(M2)) = 0;

clf;
imageplot(M1, 'Frame #1', 1,2,1);
imageplot(M2, 'Frame #2', 1,2,2);

global D;
Dt = M1-M2;
D = grad(M1);

clf;
imageplot(Dt, 'd/dt', 1,3,1);
imageplot(D(:,:,1), 'd/dx', 1,3,2);
imageplot(D(:,:,2), 'd/dy', 1,3,3);

global lambda;
lambda = .1;

tau = .2;

v = zeros(n,n,2);

U = Dt + sum(v.*D,3);

L = cat(3, div(grad(v(:,:,1))), div(grad(v(:,:,2))));

G = D.*repmat(U, [1 1 2])  - lambda * L;

v = v - tau*G;

exo1()

%% Insert your code here.

tol = 1e-5;
maxit = 200;

b = -D.*cat(3,Dt,Dt);

[v,flag,relres,it,resvec] = cgs(@callback_optical_flow,b(:),tol,maxit);
v = reshape(v, [n n 2]);

clf;
imageplot(v, '', 1,2,1);
subplot(1,2,2);
w = 12; m = ceil(n/w);
t = w/2 + ((0:m-1)*w);
[V,U] = meshgrid(t,t);
hold on;
imageplot(M1);
quiver(t,t,v(1:w:n,1:w:n,2), v(1:w:n,1:w:n,1));
axis('ij');

[Y,X] = meshgrid(1:n,1:n);
X = clamp(X+v(:,:,1),1,n);
Y = clamp(Y+v(:,:,2),1,n);

Ms = interp2( 1:n,1:n, M1, Y,X );

R0 = M2-M1;

R = Ms-M1;

v = max( [max(abs(R0(:))) max(abs(R(:)))] );
R(1)=v; R(2)=-v; R0(1)=v; R0(2)=-v;

clf;
imageplot(R0, 'Residual without flow', 1,2,1);
imageplot(R, 'Residual with flow', 1,2,2);

w = 8;

q = 4;

dq = .5;

m = ceil(n/w);

[X0,Y0,dX,dY] = ndgrid( 0:w-1, 0:w-1, -q:dq:q,-q:dq:q);
[dy,dx] = meshgrid(-q:dq:q,-q:dq:q);

F = zeros(n,n,2);

i = 3; j = 40;

x = (i-1)*w+1;
y = (j-1)*w+1;

selx = clamp( (i-1)*w+1:i*w, 1,n);
sely = clamp( (j-1)*w+1:j*w, 1,n);

X = clamp(x + X0 + dX,1,n);
Y = clamp(y + Y0 + dY,1,n);

P2 = M2(selx,sely);

P1 = interp2( 1:n,1:n, M1, Y,X );

d = sum(sum( (P1-repmat(P2,[1 1 size(P1,3) size(P1,4)])).^2 ) );

[tmp,I] = compute_min(d(:));
F(selx,sely,1) = dx(I);
F(selx,sely,2) = dy(I);

exo2()

%% Insert your code here.

clf;
imageplot(F, '', 1,2,1);
subplot(1,2,2);
t = w/2 + ((0:m-1)*w);
[V,U] = meshgrid(t,t);
hold on;
imageplot(M1);
quiver(t,t,F(1:w:n,1:w:n,2), F(1:w:n,1:w:n,1));
axis('ij');

[Y,X] = meshgrid(1:n,1:n);
X = clamp(X+F(:,:,1),1,n);
Y = clamp(Y+F(:,:,2),1,n);

Ms = interp2( 1:n,1:n, M1, Y,X );

R0 = M2-M1;

R = M2-Ms;

v = max( [max(abs(R0(:))) max(abs(R(:)))] );
R(1)=v; R(2)=-v; R0(1)=v; R0(2)=-v;

clf;
imageplot(R0, 'Residual without flow', 1,2,1);
imageplot(R, 'Residual with flow', 1,2,2);
