
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/wavelet_2_haar2d')

n = 256;
N = n*n;

name = 'hibiscus';
f = load_image(name,n);
f = rescale( sum(f,3) );

J = log2(n)-1;

haarV = @(a)[   a(1:2:length(a),:) + a(2:2:length(a),:);
                    a(1:2:length(a),:) - a(2:2:length(a),:) ]/sqrt(2);

clf;
imageplot(f,'Original image',1,2,1);
imageplot(haarV(f),'Vertical transform',1,2,2);

haarH = @(a)haarV(a')';

haar = @(a)haarH(haarV(a));

clf;
imageplot(f,'Original image',1,2,1);
subplot(1,2,2);
plot_wavelet(haar(f),log2(n)-1); title('Transformed');

j = J;

fw = f;

fw(1:2^(j+1),1:2^(j+1)) = haar(fw(1:2^(j+1),1:2^(j+1)));

exo1()

%% Insert your code here.

fprintf('Should be 0: %.3f.\n', (norm(f(:))-norm(fw(:)))/norm(f(:)));

clf;
subplot(1,2,1);
imageplot(f); title('Original');
subplot(1,2,2);
plot_wavelet(fw, 1); title('Transformed');

ihaarV = @(a,d)assign( zeros(2*size(a,1),size(a,2)), ...
        [1:2:2*size(a,1), 2:2:2*size(a,1)], [a+d; a-d]/sqrt(2), 1 );

ihaarH = @(a,d)ihaarV(a',d')';

ihaar = @(a,dH,dV,dD)ihaarV( ihaarH(a,dH), ihaarH(dV,dD) );

f1 = fw;
j = 0;

s = 1:2^j; t = 2^j+1:2^(j+1); u = 1:2^(j+1);
f1(u,u) = ihaar(f1(s,s),f1(s,t),f1(t,s),f1(t,t));

exo2()

%% Insert your code here.

fprintf('Should be 0: %.3f.\n', (norm(f-f1))/norm(f));

j = J-1;

fw1 = zeros(n); fw1(1:2^j,1:2^j) = fw(1:2^j,1:2^j);

exo3()

%% Insert your code here.

S = @(x,T) x .* (abs(x)>T);

T = .1;

fwT = S(fw,T);

clf;
subplot(1,2,1);
plot_wavelet(fw); axis('tight'); title('Original coefficients');
subplot(1,2,2);
plot_wavelet(fwT); axis('tight'); title('Thresholded coefficients');

exo4()

%% Insert your code here.
