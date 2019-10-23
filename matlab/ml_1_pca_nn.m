
addpath('toolbox_general')
addpath('solutions/ml_1_pca_nn')

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 10);

name = 'digits';
name = 'iris';
load(['ml-' name]);

A = A(randperm(size(A,1)),:);

X = A(:,1:end-1);
y = A(:,end);
y = y-min(y)+1;

[n,p] = size(X);
k = max(y);

Xm = @(X)X-repmat(mean(X,1), [size(X,1) 1]);
Cov = @(X)Xm(X)'*Xm(X);

clf;
imagesc(Cov(X));
colormap jet(256);

[U,D,V] = svd(Xm(X),'econ');

Z = Xm(X) * V;

clf; 
plot(diag(D), '.-', 'LineWidth', 2, 'MarkerSize', 30);
axis tight;
SetAR(1/2);

col = [ [1 0 0]; [0 1 0]; [0 0 1]; [0 0 0]; [0 1 1]; [1 0 1]; [1 1 0]; ...
    [1 .5 .5]; [.5 1 .5]; [.5 .5 1]  ]';
ms = 25;
clf; hold on;
lgd = {};
for i=1:min(k,size(col,2))
    I = find(y==i);
    plot(Z(I,1), Z(I,2), '.', 'Color', col(:,i), 'MarkerSize', ms);
    lgd{end+1} = num2str(i);
end
axis tight; axis equal; box on;
legend(lgd, 'Location', 'EastOutside');
SetAR(1);

clf; hold on;
for i=1:k
    I = find(y==i);
    plot3(Z(I,1), Z(I,2), Z(I,3), '.', 'Color', col(:,i), 'MarkerSize', ms);
end
view(3); axis tight; axis equal; box on;
legend(lgd, 'Location', 'EastOutside');
SetAR(1);

n0 = round(.5*n); n1 = n-n0;
X0 = X(1:n0,:);     y0 = y(1:n0);
X1 = X(n0+1:end,:); y1 = y(n0+1:end);

distmat = @(X,Z)bsxfun(@plus,dot(X',X',1)',dot(Z',Z',1))-2*(X*Z');

i = 1; x = X1(i,:);  % could be any point
D = distmat(X0,x);

[~,I] = sort(D);
ys = y(I);

R = 5; 
h = hist(ys(1:R,:), 1:k) / R;
[~,c] = max(h);
fprintf('c(x)=%d [true class=%d]\n', c, y1(i));

Rlist = round([.05 .1 .5 1]*n0); %  [5 50 100]
clf;
for i=1:length(Rlist)
    R = Rlist(i);
    h = hist(ys(1:R,:), 1:k)/R;
    subplot(length(Rlist),1,i);
    bar(1:k,h); 
    axis([0.5 k+.5 0 1]);
    set(gca, 'FontSize', 15);
end

exo1()

%% Insert your code here.

exo2()

%% Insert your code here.

if k>=4
ksvg = k; Xsvg = X; ysvg = y;
k = 3;
I = find(y<=k);
X = X(I,:); y = y(I);
n = length(I);
end

[U,D,V] = svd(Xm(X),'econ');
Z = Xm(X) * V;

I = randperm(n); I = I(1:k);
C = X(I,:);

D = distmat(X,C);
[~,yb] = min(D, [], 2);

clf;
hold on;
for i=1:k
    I = find(yb==i);
    plot(Z(I,1), Z(I,2), '.', 'Color', col(:,i), 'MarkerSize', 25);
end
CV = (C-repmat(mean(X,1), [k 1]))*V;
for i=1:k
    plot(CV(i,1), CV(i,2), 'o', 'MarkerFaceColor', col(:,i), 'MarkerSize', 12, 'MarkerEdgeColor', 'k');
end
axis tight; axis equal; axis off;
SetAR(1);

for l=1:k
    C(l,:) = mean( X(yb==l,:), 1 );
end

exo3()

%% Insert your code here.

clf
for l=1:k
    I = find(yb==l);
    h = hist(y(I),1:k); h = h/sum(h);
    subplot(k,1,l);
    bar(1:k,h); 
    axis([0.5 k+.5 0 1]);
    set(gca, 'FontSize', 10);
end

exo4()

%% Insert your code here.
<script>
  $(document).ready(function(){
      $('div.prompt').hide();
  });
</script>