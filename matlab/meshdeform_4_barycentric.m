
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/meshdeform_4_barycentric')

cage = 'V';

delta = .03;
rho = .2;
eta = .05;
x1 = .5-rho-eta/2; x2 = .5-eta/2;
x3 = .5+eta/2; x4 = .5+rho+eta/2;
switch cage
    case 'L'
        V = [[delta;delta] [1-delta;delta] [1-delta;delta+rho] [delta+rho;delta+rho] [delta+rho;1-delta] [delta;1-delta]];
    case 'U'
        V = [[x1;delta] [x4;delta] [x4;1-delta] [x3;1-delta] ...
            [x3;delta+rho] [x2;delta+rho] [x2;1-delta] [x1;1-delta] ];
    case 'V'
        V = [[x1;delta] [x4;delta] [x4;1-delta] [x3;1-delta] ...
            [.5;delta+rho] [x2;1-delta] [x1;1-delta] ];
end
k = size(V,2);

n = 200;
x = linspace(0,1,n);
[Y,X] = meshgrid(x,x);

S = 1 - inpolygon(X,Y,V(1,:),V(2,:));

m  = 30;
[XT,YT] = meshgrid((0:n-1)/n,(0:n-1)/n);
T = mod( floor(XT*m)+floor(YT*m),2 );
T(S==1) = 1;

lw = 3; ms = 25;
clf; hold on;
plot_surf_texture(cat(3,X,Y,zeros(n)), T');
h = plot(V(1,[1:end 1]), V(2,[1:end 1]), 'r.-');
set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms);
view(2); axis('off'); axis('equal');

dotp = @(a,b)sum(a.*b);
crossp = @(a,b)a(1,:).*b(2,:)-a(2,:).*b(1,:);
normalize = @(a)a./repmat(sqrt(sum(a.^2)), [2 1]);

W = [X(:)';Y(:)'];

C = zeros(n,n,k);

i = 4;
vi = V(:,i);

U = repmat(vi,[1 n^2])-W;
nb = normalize( U );

d = sqrt( sum(U.^2) );
for j=mod([i-2,i],k)+1
    % point
    vj = V(:,j);
    na = normalize( repmat(vj,[1 n^2])-W );
    % angle
    dp = dotp(na,nb);
    theta = acos(clamp(dp,-1,1));
    % add tangent of half angle
    C(:,:,i) = C(:,:,i) + reshape( tan(theta/2) ./ d, [n n]);
end

exo1()

%% Insert your code here.

C = C ./ repmat( sum(C,3), [1 1 k] );

i = 4;
c = abs(C(:,:,i))+1e-3;
c(S==1) = 0;

nl = 15;
t = linspace(0,1,n);
B = display_shape_function(c');
clf; hold on;
imagesc(t,t,B); axis('image'); axis('off');
contour(t,t,c',nl, 'k');
colormap jet(256);
h = plot(V(1,[1:end 1]), V(2,[1:end 1]), 'r.-');
set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms);

applyinterp = @(C,x)sum(repmat(reshape(x(:), [1 1 k]),[n n 1]).*C,3);

clf;
imageplot(applyinterp(C,V(1,:)), 'Should be X', 1,2,1);
imageplot(applyinterp(C,V(2,:)), 'Should be Y', 1,2,2);
colormap jet(256);

V2 = V;
switch cage,
    case 'L'
        V2(:,4) = [1-delta;1-delta];
    case {'U' 'V'}
        V2(:,3) = [1-delta;delta];
        V2(:,4) = [1-delta;delta+rho];
end

rho = .7;
V1 = V*(1-rho) + V2*rho;

X1 = applyinterp(C,V1(1,:));
Y1 = applyinterp(C,V1(2,:));

clf; hold on;
plot_surf_texture(cat(3,X1,Y1,zeros(n)), T');
h = plot(V1(1,[1:end 1]), V1(2,[1:end 1]), 'r.-');
set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms);
view(2); axis('off'); axis('equal');
axis([0,1,0,1]);

Vi = round(V*n);
bound = Vi(:,1);
loc = 1;
abscur = 0;
for i=1:k
    j = mod(i,k)+1;
    d = round(norm(Vi(:,i)-Vi(:,j)));
    t = repmat((1:d)/d, [2 1]);
    bound = [bound repmat(Vi(:,i),[1 d]).*(1-t) + repmat(Vi(:,j),[1 d]).*t];
    abscur = [abscur abscur(end)+(1:d)/d];
    loc(end+1) = loc(end)+d;
end
q = size(bound,2);
bound = round(bound);
I = bound(1,:) + (bound(2,:)-1)*n;

C = zeros(n,n,k);

i = 4;
u = zeros(k+1,1);
u(i) = 1;
if i==1
    u(end)=1;
end
u = interp1(0:k,u, abscur);

Ci = zeros(n);

sel1 = [2:n 1];
sel2 = [n 1:n-1];
Ci = ( Ci(sel1,:) + Ci(:,sel1) + Ci(sel2,:) + Ci(:,sel2) )/4;

Ci(I) = u;

clf;
imageplot(Ci');
colormap jet(256);

exo2()

%% Insert your code here.

c = C(:,:,i); c(S==1) = 0;
nl = 15;
t = linspace(0,1,n);
B = display_shape_function(c');
clf; hold on;
imagesc(t,t,B); axis('image'); axis('off');
contour(t,t,c',nl, 'k');
colormap jet(256);
h = plot(V(1,[1:end 1]), V(2,[1:end 1]), 'r.-');
set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms);

exo3()

%% Insert your code here.

clf;
A = applyinterp(C,V(1,:)); A(S==1) = 0;
imageplot(A, 'Should be X', 1,2,1);
A = applyinterp(C,V(2,:)); A(S==1) = 0;
imageplot(A, 'Should be Y', 1,2,2);
colormap jet(256);

rho = .7;
V1 = V*(1-rho) + V2*rho;

X1 = applyinterp(C,V1(1,:));
Y1 = applyinterp(C,V1(2,:));
X1(S==1) = Inf; Y1(S==1) = Inf;

clf; hold on;
plot_surf_texture(cat(3,X1,Y1,zeros(n)), T');
h = plot(V1(1,[1:end 1]), V1(2,[1:end 1]), 'r.-');
set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms);
view(2); axis('off'); axis('equal');
axis([0,1,0,1]);

N = V(:,[2:end 1]) - V;
N = N ./ repmat( sqrt(sum(N.^2)), [2 1] );
N = -[-N(2,:); N(1,:)];

clf; hold on;
h = plot(V(1,[1:end 1]), V(2,[1:end 1]), 'r.-');
set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms);
rho = .1;
for i=1:k
    j = mod(i,k)+1;
    a = mean( V(:,[i,j]), 2);
    h = plot( [a(1) a(1)+rho*N(1,i)], [a(2) a(2)+rho*N(2,i)], 'k' );
    set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms);
end
axis square; axis off;

C = zeros(n,n,k);
D = zeros(n,n,k);

i = 4;
j = mod(i,k)+1;
vi = V(:,i);
vj = V(:,j);
ni = -N(:,i);

a = repmat(vj - vi,[1 n^2]);
b = repmat(vi,[1 n^2]) - W;
Q = sum(a.^2); s = sum(b.^2); R = 2*sum(a.*b);
na = sqrt(sum(a.^2));
BA = na .* sum( b .* repmat(ni,[1 n^2]) );
SRT = sqrt( 4*s.*Q - R.^2 );
L0 = log(s); L1 = log(s+Q+R);
A0 = atan(      R ./SRT ) ./ SRT;
A1 = atan( (2*Q+R)./SRT ) ./ SRT;
A10 = A1 - A0;
L10 = L1 - L0;

d = - na .* ( (4*s-(R.^2)./Q) .* A10 + R./(2*Q).*L10 + L1 - 2 )  / (4*pi) ;
cj = - BA .* ( L10./(2*Q) - A10 .* (  R./Q) ) / (2*pi);
ci = + BA .* ( L10./(2*Q) - A10 .* (2+R./Q) ) / (2*pi);
D(:,:,i) = reshape(d,n,n);
C(:,:,i) = C(:,:,i) + reshape(ci,n,n);
C(:,:,j) = C(:,:,j) + reshape(cj,n,n);

exo4()

%% Insert your code here.

i = 4;
c = rescale(C(:,:,i)); c(S==1) = 0;
nl = 15;
t = linspace(0,1,n);
B = display_shape_function(c');
clf; hold on;
imagesc(t,t,B); axis('image'); axis('off');
contour(t,t,c',nl, 'k');
colormap jet(256);
h = plot(V(1,[1:end 1]), V(2,[1:end 1]), 'r.-');
set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms);

c = D(:,:,i); c(S==1) = 0;
nl = 15;
t = linspace(0,1,n);
B = display_shape_function(c');
clf; hold on;
imagesc(t,t,B); axis('image'); axis('off');
contour(t,t,c',nl, 'k');
colormap jet(256);
h = plot(V(1,[1:end 1]), V(2,[1:end 1]), 'r.-');
set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms);

x = applyinterp(C,V(1,:)) + applyinterp(D,N(1,:));
y = applyinterp(C,V(2,:)) + applyinterp(D,N(2,:));
x(S==1) = 0; y(S==1) = 0;
clf;
imageplot(x, 'Should be X', 1,2,1);
imageplot(y, 'Should be Y', 1,2,2);
colormap jet(256);

rho = .7;
V1 = V*(1-rho) + V2*rho;

N1 = V1(:,[2:end 1]) - V1;
N1 = N1 ./ repmat( sqrt(sum(N1.^2)), [2 1] );
N1 = -[-N1(2,:); N1(1,:)];

s = sqrt( sum( (V1(:,[2:end 1]) - V1).^2 ) ./ sum( (V(:,[2:end 1]) - V).^2 ) );
s = repmat(reshape(s, [1 1 k]), [n n 1]);

X1 = applyinterp(C,V1(1,:)) + applyinterp(s.*D,N1(1,:));
Y1 = applyinterp(C,V1(2,:)) + applyinterp(s.*D,N1(2,:));
X1(S==1) = Inf; Y1(S==1) = Inf;

clf; hold on;
plot_surf_texture(cat(3,X1,Y1,zeros(n)), T');
h = plot(V1(1,[1:end 1]), V1(2,[1:end 1]), 'r.-');
set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms);
view(2); axis('off'); axis('equal');
axis([0,1,0,1]);

exo5()

%% Insert your code here.

exo6()

%% Insert your code here.
