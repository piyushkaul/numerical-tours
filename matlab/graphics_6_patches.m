
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/graphics_6_patches')

w = 5;

q = 1000;

name = 'corral';
n = 128;
M0 = load_image(name,n);
M0 = rescale( crop(M0,n) );

p = n-w+1;

[Y,X] = meshgrid(1:p,1:p);

[dY,dX] = meshgrid(0:w-1,0:w-1);

X = reshape(X, [1 1 p p]);
Y = reshape(Y, [1 1 p p]);
X = repmat(X, [w w 1 1]) + repmat(dX, [1 1 p p]);
Y = repmat(Y, [w w 1 1]) + repmat(dY, [1 1 p p]);

P0 = M0(X + (Y-1)*n);
P0 = reshape(P0,w,w,p*p);

sel = randperm(size(P0,3)); sel = sel(1:q);
P0 = P0(:,:,sel);

n = 128;
M = rand(n);

ofx = 2;
ofy = 1;

[Y,X] = meshgrid(1:w:n, 1:w:n);
p = size(X,1);
[dY,dX] = meshgrid(0:w-1,0:w-1);
X = reshape(X, [1 1 p p]);
Y = reshape(Y, [1 1 p p]);
X = repmat(X, [w w 1 1]) + repmat(dX, [1 1 p p]);
Y = repmat(Y, [w w 1 1]) + repmat(dY, [1 1 p p]);

Xs = mod(X+ofx-1, n)+1;
Ys = mod(Y+ofy-1, n)+1;

P = M(Xs + (Ys-1)*n);

for i=1:p*p
    % distance to current patch
    d = sum(sum( (P0 - repmat(P(:,:,i), [1 1 q])).^2 ) );
    % best match
    [tmp,s] = min(d);
    % replace the patch
    P(:,:,i) = P0(:,:,s);
end

Mp = M;
Mp(Xs + (Ys-1)*n) = P;

clf;
imageplot(M,'Input', 1,2,1);
imageplot(Mp,'Projected', 1,2,2);

noffs = 10;

sel = randperm(w*w); sel = sel(1:noffs);
OffX = dX(sel); OffY = dY(sel);

M1 = zeros(n);
for j=1:noffs
    ofx = OffX(j);
    ofy = OffY(j);
    % shift locations
    Xs = mod(X+ofx-1, n)+1;
    Ys = mod(Y+ofy-1, n)+1;
    % extract patch
    P = M(Xs + (Ys-1)*n);
    % replace by closest patch
    for i=1:p*p
        d = sum(sum( (P0 - repmat(P(:,:,i), [1 1 q])).^2 ) );
        [tmp,s] = min(d);
        P(:,:,i) = P0(:,:,s);
    end
    % reconstruct the image.
    M1(Xs + (Ys-1)*n) = M1(Xs + (Ys-1)*n) + P;
end
M1 = M1 / noffs;

M1 = perform_hist_eq(M1,M0);

clf;
imageplot(M,'Input', 1,2,1);
imageplot(M1,'Projected', 1,2,2);

exo1()

%% Insert your code here.

clf;
imageplot(M0,'Exemplar', 1,2,1);
imageplot(M,'Synthesized', 1,2,2);

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.

name = 'hair';
n = 128;
Ma = load_image(name);
Ma = rescale( crop(Ma,n) );

h = 9;

nh = 20;
mask = zeros(n);
[V,U] = meshgrid(1:n,1:n);
for i=1:nh
    % location of the hole
    x = floor(rand(2,1)*(n-1))+1;
    d = (U-x(1)).^2 + (V-x(2)).^2;
    mask( d<=h^2 ) = 1;
end

M0 = Ma;
M0(mask==1) = 0;

clf;
imageplot(Ma, 'Original', 1,2,1);
imageplot(M0, 'To inpaint', 1,2,2);

p = n-w+1;
[Y,X] = meshgrid(1:p,1:p);
[dY,dX] = meshgrid(0:w-1,0:w-1);
X = reshape(X, [1 1 p p]);
Y = reshape(Y, [1 1 p p]);
X = repmat(X, [w w 1 1]) + repmat(dX, [1 1 p p]);
Y = repmat(Y, [w w 1 1]) + repmat(dY, [1 1 p p]);
P0 = M0(X + (Y-1)*n);
P0 = reshape(P0,w,w,p*p);

I = find( min(min(P0,[],1),[],2)~=0 );
P0 = P0(:,:,I);

q = 1000;

sel = randperm(size(P0,3)); sel = sel(1:q);
P0 = P0(:,:,sel);

M = M0;
I = find(mask==1);
M(I) = rand(length(I),1);

ofx = 2; ofy = 1;
[Y,X] = meshgrid(1:w:n, 1:w:n);
p = size(X,1);
[dY,dX] = meshgrid(0:w-1,0:w-1);
X = reshape(X, [1 1 p p]);
Y = reshape(Y, [1 1 p p]);
X = repmat(X, [w w 1 1]) + repmat(dX, [1 1 p p]);
Y = repmat(Y, [w w 1 1]) + repmat(dY, [1 1 p p]);
Xs = mod(X+ofx-1, n)+1;
Ys = mod(Y+ofy-1, n)+1;
P = M(Xs + (Ys-1)*n);
Pmask = M(Xs + (Ys-1)*n);

for i=1:p*p
    if sum(sum(Pmask(:,:,i)))>0
        % project only a patch crossing the hole.
        % distance to current patch
        d = sum(sum( (P0 - repmat(P(:,:,i), [1 1 q])).^2 ) );
        % best match
        [tmp,s] = min(d);
        % replace the patch
        P(:,:,i) = P0(:,:,s);
    end
end

Mp = M;
Mp(Xs + (Ys-1)*n) = P;

Mp(mask==0) = M0(mask==0);

clf;
imageplot(M,'Input', 1,2,1);
imageplot(Mp,'Projected', 1,2,2);

exo6()

%% Insert your code here.

clf;
imageplot(M0, 'To inpaint', 1,2,1);
imageplot(M, 'Inpainted', 1,2,2);

exo7()

%% Insert your code here.
