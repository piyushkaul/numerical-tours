
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/multidim_1_color')

n = 256;
N = n*n;

name = 'hibiscus';
f = rescale( load_image(name,n) );

R = cat(3, f(:,:,1), zeros(n), zeros(n));
G = cat(3, zeros(n), f(:,:,2), zeros(n));
B = cat(3, zeros(n), zeros(n), f(:,:,3));
clf;
imageplot({f R G B}, ...
        { 'f' 'R (Red)' 'G (green)' 'B (blue)'}, 2, 2);

clf;
imageplot({f mean(f,3)}, {'f' 'L'});

f1 = cat(3, f(:,:,1),     f(:,:,2)*0+1, f(:,:,3)*0+1);
f2 = cat(3, f(:,:,1)*0+1, f(:,:,2)    , f(:,:,3)*0+1);
f3 = cat(3, f(:,:,1)*0+1, f(:,:,2)*0+1, f(:,:,3));
clf;
imageplot({f f1 f2 f3}, ...
        { 'f' 'C' 'f' 'Y'}, 2, 2);

T = [.299 .587 .114; ...
    -.14713 -.28886 .436; ...
    .615 -.51499 -.10001]';

applymat = @(f,T)reshape( reshape(f, [n*n 3])*T, [n n 3] );
rgb2yuv  = @(f)applymat(f,T);

U = rgb2yuv(f);
clf;
imageplot(U(:,:,1), 'Y', 1,3,1);
imageplot(U(:,:,2), 'U', 1,3,2);
imageplot(U(:,:,3), 'V', 1,3,3);

U1 = U;
U1(:,:,2:3) = U1(:,:,2:3)/2;

exo1()

%% Insert your code here.

Value = @(f)sum(f, 3) / sqrt(3);

A = @(f)( f(:,:,2)-f(:,:,3) )/sqrt(2);
B = @(f)( 2*f(:,:,1) - f(:,:,2) - f(:,:,3) )/sqrt(6);

T = [   1/sqrt(3) 1/sqrt(3) 1/sqrt(3); ...
        0 1/sqrt(2) -1/sqrt(2); ...
        2/sqrt(6) -1/sqrt(6) -1/sqrt(6)];

Saturation = @(f)sqrt( A(f).^2 + B(f).^2 );
Hue = @(f)atan2(B(f),A(f));

rgb2hsv1 = @(f)cat(3, Hue(f), Saturation(f), Value(f));

g = rgb2hsv1(f);

clf;
imageplot({g(:,:,1) g(:,:,2) g(:,:,3)}, {'H' 'S' 'V'}, 1,3);

exo2()

%% Insert your code here.

m = mean(mean(f,1), 2);

X = reshape( f - repmat(m, [n n 1]), [n*n 3] );

C = (X'*X)/N;

[V,D] = eig(C); D = diag(D);
[D,I] = sort(D, 'descend'); V = V(:,I);

rgb2pca = @(f,V,m)applymat(f - repmat(m, [n n 1]),V);
g = rgb2pca(f,V,m);

clf;
imageplot({g(:,:,1) g(:,:,2) g(:,:,3)}, {'PCA_1' 'PCA_2' 'PCA_3'}, 1,3);

g1 = g;
g1(:,:,2:3) = g1(:,:,2:3)/2;

exo3()

%% Insert your code here.

channel = @(f,i)reshape(f(:,:,i), [N 1]);

Q = 60;

clf; c = {'r' 'g' 'b'}; lgd = {'R' 'G' 'B'};
for i=1:3
    subplot(3,1,i); 
    [h,t] = hist(channel(f,i), Q);
    bar(t,h*Q/N, c{i}); axis('tight');
    legend(lgd{i});
end

g = rgb2hsv1(f);
clf; c = {'k' 'k' 'k'}; lgd = {'H' 'S' 'V'};
for i=1:3
    subplot(3,1,i); 
    [h,t] = hist(channel(g,i), Q);
    bar(t,h*Q/N, c{i}); axis('tight');
    legend(lgd{i});
end

P = 5000;

H = reshape(f, [n*n 3]);
sel = randperm(n*n); sel = sel(1:P);
H = H(sel,:);

plotp = @(x,col)plot3(x(1,:)', x(2,:)', x(3,:)', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', 2);
clf; hold on;
for i=1:P
    plotp(H(i,:)', H(i,:));
end
view(3);

quantize = @(A,Q)1+round((Q-1)*A);
J = @(I,Q)I(1,:)' + Q*(I(2,:)'-1);
hist2d = @(f,Q)reshape( accumarray(J(quantize(f,Q),Q), ones(1,N), [Q*Q 1], @sum), [Q Q]);

Q = 60;

func = @(a)log(a+3);
X = reshape(f(:,:,1:2), [n*n 2])';
clf;
imageplot( func(hist2d(X,Q)) );

sigma = .13;

f1 = f + randn(n,n,3)*sigma;

clf;
imageplot(f, 'f', 1,2,1);
imageplot(clamp(f1), 'f_1', 1,2,2);

exo4()

%% Insert your code here.

exo5()

%% Insert your code here.

exo6()

%% Insert your code here.
