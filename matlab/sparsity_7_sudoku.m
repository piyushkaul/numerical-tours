
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/sparsity_7_sudoku')

n = 9;

x = floor(rand(n)*n)+1;

U = repmat( reshape(1:n, [1 1 n]), n,n );

encode = @(x)double( repmat(x, [1 1 n])==U );
X = encode(x);

[tmp,x1] = min( abs(X-1), [], 3 );

disp(['Should be 0: ' num2str(norm(x-x1,'fro')) '.']);

i = 4; j = 4;
Z = zeros(n,n,n);
Z(i,j,:) = 1;

Aenc = [];
Aenc(end+1,:) = Z(:)';

disp(['Should be 1: ' num2str(Aenc*X(:)) '.']);

exo1()

%% Insert your code here.

disp(['Should be 0: ' num2str(norm(Aenc*X(:)-1)) '.']);

x = [8 1 9 6 7 4 3 2 5; 
     5 6 3 2 8 1 9 4 7; 
     7 4 2 5 9 3 6 8 1; 
     6 3 8 9 4 5 1 7 2; 
     9 7 1 3 2 8 4 5 6; 
     2 5 4 1 6 7 8 9 3;
     1 8 5 7 3 9 2 6 4;
     3 9 6 4 5 2 7 1 8; 
     4 2 7 8 1 6 5 3 9]

X = encode(x);

i=3; k=5;
Z = zeros(n,n,n);
Z(i,:,k) = 1;

Arow = [];
Arow(end+1,:) = Z(:)';

disp(['Should be 1: ' num2str(Arow*X(:)) '.']);

exo2()

%% Insert your code here.

disp(['Should be 0: ' num2str(norm(Arow*X(:)-1)) '.']);

exo3()

%% Insert your code here.

disp(['Should be 0: ' num2str(norm(Acol*X(:)-1)) '.']);

p = sqrt(n);

k = 1;
Z = zeros(n,n,n);
Z(1:p,1:p,k) = 1;

Ablock = [];
Ablock(end+1,:) = Z(:)';

disp(['Should be 1: ' num2str(Ablock*X(:)) '.']);

exo4()

%% Insert your code here.

disp(['Should be 0: ' num2str(norm(Ablock*X(:)-1)) '.']);

x1 = [0 1 0 0 0 0 3 0 0; 
     0 0 3 0 8 0 0 4 0; 
     7 0 2 0 0 3 0 0 1; 
     0 3 0 9 4 0 1 0 0; 
     9 0 0 0 0 0 0 0 6; 
     0 0 4 0 6 7 0 9 0;
     1 0 0 7 0 0 2 0 4;
     0 9 0 0 5 0 7 0 0; 
     0 0 7 0 0 0 0 3 0];

[I,J] = ind2sub( [n n], find(x1(:)~=0) );
v = x1(x1(:)~=0);

i = 1;
Z = zeros(n,n,n);
Z(I(i), J(i), v(i)) = 1;

Ainp = [];
Ainp(end+1,:) = Z(:)';

exo5()

%% Insert your code here.

X1 = encode(x1);
disp(['Should be 0: ' num2str(norm(Ainp*X1(:)-1)) '.']);

A = [Aenc; Arow; Acol; Ablock; Ainp];

pA = pinv(A);

exo6()

%% Insert your code here.

projector = @(u)reshape( u(:) - pA*(A*u(:)-1), [n n n]);

Xproj = projector( zeros(n,n,n) );

d = projector(Xproj)-Xproj;
disp(['Should be 0: ' num2str(norm(d(:), 'fro')) '.']);

clf;
hist(Xproj(:), 30);
axis('tight');

[tmp,xproj] = min( abs(Xproj-1), [], 3 );
Xproj1 = encode(xproj); 
disp(['Number of violated constraints: ' num2str(sum(A*Xproj1(:)~=1)) '.']);

Xproj = zeros(n,n,n);

Xproj = max( projector(Xproj),0 );

exo7()

%% Insert your code here.

clf;
hist(Xproj(:), 30);
axis('tight');

[tmp,xproj] = min( abs(Xproj-1), [], 3 );
Xproj1 = encode(xproj); 
disp(['Number of violated constraints: ' num2str(sum(A*Xproj1(:)~=1)) '.']);

exo8()

%% Insert your code here.

x1 = [0 0 3 0 0 9 0 8 1;
     0 0 0 2 0 0 0 6 0;
     5 0 0 0 1 0 7 0 0;
     8 9 0 0 0 0 0 0 0;
     0 0 5 6 0 1 2 0 0;
     0 0 0 0 0 0 0 3 7;
     0 0 9 0 2 0 0 0 8;
     0 7 0 0 0 4 0 0 0;
     2 5 0 8 0 0 6 0 0];

exo9()

%% Insert your code here.

Xproj = projector( zeros(n,n,n) );
disp(['Sparsity: ' num2str(sum(Xproj(:)~=0)) ' (optimal: ' num2str(n*n) ').']);

disp(['L1 norm: ' num2str(norm(Xproj(:),1)) ' (optimal: ' num2str(n*n) ').']);

solvel1 = @(A)reshape(perform_solve_bp(A, ones(size(A,1),1), n^3, 30, 0, 1e-10), [n n n]);

Xbp = solvel1(A);

disp(['L1 norm: ' num2str(norm(Xbp(:),1)) ' (optimal: ' num2str(n*n) ').']);

[tmp,xbp] = min( abs(Xbp-1), [], 3 );
Xbp1 = encode(xbp); 
disp(['Number of violated constraints: ' num2str(sum(A*Xbp1(:)~=1)) '.']);

alpha = 0;

epsilon = 0.1;

u = ones(n,n,n);

Xrw = solvel1( A*diag(u(:)) ) .* u;

u = (abs(Xrw).^(1-alpha)+epsilon);

exo10()

%% Insert your code here.

hist(Xrw(:), 30);

x1 = [1 0 0 0 0 7 0 9 0;
      0 3 0 0 2 0 0 0 8;
      0 0 9 6 0 0 5 0 0; 
      0 0 5 3 0 0 9 0 0; 
      0 1 0 0 8 0 0 0 2;
      6 0 0 0 0 4 0 0 0; 
      3 0 0 0 0 0 0 1 0; 
      0 4 0 0 0 0 0 0 7;
      0 0 7 0 0 0 3 0 0];

exo11()

%% Insert your code here.

exo12()

%% Insert your code here.

exo13()

%% Insert your code here.

exo14()

%% Insert your code here.
