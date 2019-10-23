
p = 500; % size of the image for sampling, m=p*p
t = linspace(0,1,p);
[V,U] = meshgrid(t,t);
Y = [U(:)';V(:)'];

n = 30;
X = .5+.5i + exp(1i*pi/4) * 1*( .1*(rand(1,n)-.5)+1i*(rand(1,n)-.5) );
X = [real(X);imag(X)];
a = ones(n,1)/n;

Z = {[.6 .9] [.4 .1]};  % mean
S = [.07 .09]; % variance
W = [.5 .5]; 
b = zeros(p); % ones(p)*.01;
for i=1:length(Z)
    z = Z{i};
    s = S(i);
    b = b  + W(i) * exp( (-(U-z(1)).^2-(V-z(2)).^2)/(2*s^2) );
end
b = b/sum(b(:));

Col = rand(n,3);
clf; hold on;
imagesc(t,t,-b);
s = 60*ones(n,1); % size
scatter( X(2,:), X(1,:), s, .8*Col, 'filled' );
axis equal; axis([0 1 0 1]); axis off;
colormap gray(256);

f = zeros(n,1);

D = sum(Y.^2)' + sum(X.^2) - 2*Y'*X - f(:)';
[fC,I] = min(D,[], 2);
I = reshape(I, [p p]);

OT = sum(f.*a) + sum(fC.*b(:));

clf;
hold on;
imagesc(t,t,I); axis image; axis off;
contour(t,t,I, -.5:n+.5, 'k', 'LineWidth', 2);
colormap(Col);
scatter( X(2,:), X(1,:), (1+(f-min(f))*10)*100, Col*.8, 'filled' );

r = sum(sum( ( I==reshape(1:n,[1 1 n]) ) .* b, 1 ),2); r = r(:);
nablaE = a-r(:);

exo1()

%% Insert your code here.

clf;
plot(1:niter, E, 'LineWidth', 2);
axis tight;

f = 0;

i = (rand<W(1))+1; % select one of the two Gaussian
y = [S(i) * randn + Z{i}(1);S(i) * randn + Z{i}(2)];

[~,i] = min( sum(y.^2)' + sum(X.^2) - 2*y'*X - f(:)' );

nablaEy = a; nablaEy(i) = nablaEy(i) - 1;

exo2()

%% Insert your code here.

clf;
plot(1:niter, E, 'LineWidth', 2);
axis tight;

X1 = X;

D = sum(Y.^2)' + sum(X1.^2) - 2*Y'*X1;
[fC,I] = min(D,[], 2);
I = reshape(I, [p p]);

A = ( I==reshape(1:n,[1 1 n]) ) .* b;
B = ( I==reshape(1:n,[1 1 n]) ) .* b .* ( U+1i*V );
X1 = sum(sum(B,1),2) ./ sum(sum(A,1),2);
X1 = [real(X1(:))';imag(X1(:))'];

exo3()

%% Insert your code here.

clf;
plot(1:niter, E, 'LineWidth', 2);
axis tight;
