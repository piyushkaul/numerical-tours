
addpath('toolbox_general')
addpath('solutions/ml_3_classification')

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 10);
Xm = @(X)X-repmat(mean(X,1), [size(X,1) 1]);
Cov = @(X)Xm(X)'*Xm(X);
dotp = @(u,v)sum(u(:).*v(:));

n = 1000; % number of sample
p = 2; % dimensionality
omega = [1 .5]*5; % offset 
X = [randn(n/2,2); randn(n/2,2)+ones(n/2,1)*omega];
y = [ones(n/2,1);-ones(n/2,1)];

options.disp_dim = 2;
options.ms = 5;
clf; plot_multiclasses(X,y,options);
SetAR(1);

t = linspace(-3,3,255)';
clf;
plot(t, [t>0, log(1+exp(t)), max(t,0)], 'LineWidth', 2 );
axis tight;
legend('Binary', 'Logistic', 'Hinge', 'Location', 'NorthWest');
SetAR(1/2);

L = @(s,y)1/n * sum( log( 1 + exp(-s.*y) ) );
E = @(w,X,y)L(X*w,y);

theta = @(v)1 ./ (1+exp(-v));
nablaL = @(s,r)- 1/n * y.* theta(-s.*y);
nablaE = @(w,X,y)X'*nablaL(X*w,y);

AddBias = @(X)[X ones(size(X,1),1)];

w = zeros(p+1,1);

tau = .8; % here we are using a fixed tau
w = w - tau * nablaE(w,AddBias(X),y);

exo1()

%% Insert your code here.

q = 201; tx = linspace(min(X(:,1)),max(X(:,1)),q); ty = linspace(min(X(:,2)),max(X(:,2)),q);
[B,A] = meshgrid( ty,tx );
G = [A(:), B(:)];

Theta = theta(AddBias(G)*w);
Theta = reshape(Theta, [q q]);

clf; hold on;
imagesc(tx,ty, Theta');
options.disp_legend = 1;
plot_multiclasses(X,y,options);
SetAR(1);

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

distmat = @(X,Z)bsxfun(@plus,dot(X',X',1)',dot(Z',Z',1))-2*(X*Z');

kappa = @(X,Z,sigma)exp( -distmat(X,Z)/(2*sigma^2) );

n = 1000; p = 2;
t = 2*pi*rand(n/2,1);
R = 2.5; 
r = R*(1 + .2*rand(n/2,1)); % radius
X1 = [cos(t).*r, sin(t).*r];
X = [randn(n/2,2); X1];
y = [ones(n/2,1);-ones(n/2,1)];

options.disp_dim = 2;
options.disp_legend = 1;
clf; plot_multiclasses(X,y,options);
axis off;
SetAR(1);

sigma = 1;
K = kappa(X,X,sigma);

F = @(h,K,y)L(K*h,y);
nablaF = @(h,K,y)K'*nablaL(K*h,y);

exo4()

%% Insert your code here.

q = 201;
tmax = 3.5;
t = linspace(-tmax,tmax,q);
[B,A] = meshgrid( t,t );
G = [A(:), B(:)];
Theta = reshape( theta(kappa(G,X,sigma)*h) , [q,q]);

clf; hold on;
imagesc(t,t, Theta');
options.disp_legend = 0;
plot_multiclasses(X,y,options);
colormap jet(256); caxis([0 1]);
axis off;

exo5()

%% Insert your code here.

exo6()

%% Insert your code here.

LSE = @(S)log( sum(exp(S), 2) );

max2 = @(S)repmat(max(S,[],2), [1 size(S,2)]);
LSE = @(S)LSE( S-max2(S) ) + max(S,[],2);

SM = @(S)exp(S) ./ repmat( sum(exp(S),2), [1 size(S,2)]);

SM = @(S)SM(S-max2(S));

name = 'digits';
load(['ml-' name]);
A = A(randperm(size(A,1)),:);
X = A(:,1:end-1); y = A(:,end);

[n,p] = size(X);
CL = unique(y); % list of classes.
k = length(CL);

q = 5;
clf;
for i=1:k
    I = find(y==CL(i));
    for j=1:q
        f = reshape(X(I(j),:), sqrt(p)*[1 1])';
        subplot(q,k, (j-1)*k+i );
        imagesc(-f); axis image; axis off;
    end
end
colormap gray(256);

options.disp_dim = 2;
options.ms = 5;
options.disp_legend = 1;
clf; plot_multiclasses(X,y,options); SetAR(1);

options.disp_dim = 3;
clf; plot_multiclasses(X,y,options); SetAR(1);

D = double( repmat(CL(:)', [n,1]) == repmat(y, [1,k]) );

E = @(W)1/n*( sum(LSE(X*W)) - dotp(X*W,D)  );

nablaE = @(W)1/n * X'* ( SM(X*W) -  D );

exo7()

%% Insert your code here.

[U,D,V] = svd(Xm(X),'econ');
Z = Xm(X) * V;
M = max(abs(Z(:)));
q = 201;
t = linspace(-M,M,q);
[B,A] = meshgrid(t,t);
G = zeros(q*q,p);
G(:,1:2) = [A(:), B(:)];
G = G*V' + repmat( mean(X,1), [q*q 1] );

Theta = SM(G*W);
Theta = reshape(Theta, [q q k]);

clf;
for i=1:k
    subplot(3,4,i);
    imagesc(Theta(:,:,i)');
    title(['Class ' num2str(i)]);
    axis image; axis off;
    colormap jet(256); SetAR(1);
end

col = [ [1 0 0]; [0 1 0]; [0 0 1]; [0 0 0]; [0 1 1]; [1 0 1]; [1 1 0]; ...
    [1 .5 .5]; [.5 1 .5]; [.5 .5 1]  ]';
R = zeros(q,q,3);
for i=1:k
    for a=1:3
        R(:,:,a) = R(:,:,a) + Theta(:,:,i) .* col(a,i);
    end
end

options.disp_dim = 2;
options.ms = 3;
options.disp_legend = 1;
clf; hold on;
imagesc(t, t, permute(R, [2 1 3])); 
plot_multiclasses(X,y,options);
axis off; SetAR(1);

exo8()

%% Insert your code here.
<script>
  $(document).ready(function(){
      $('div.prompt').hide();
  });
</script>