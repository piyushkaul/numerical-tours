
warning off
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/coding_2_entropic')
warning on

p = 0.1;

n = 512;

x = (rand(n,1)>p)+1;

h = hist(x, [1 2]);
h = h/sum(h);
disp(strcat(['Empirical p=' num2str(h(1)) '.']));

e = - sum( h .* log2( max(h,1e-20) ) );
disp( strcat(['Entropy=' num2str(e)]) );

h = [.1 .15 .4 .15 .2];

m = length(h);
T = cell(0); % create an empty cell
for i=1:m
    T = cell_set(T,i,i);
end

p = h;

while length(p)>1   
    % sort in decaying order the probabilities
    [v,I] = sort(p);
    if v(1)>v(length(v))
        v = reverse(v); I = reverse(I);
    end 
    q = sum(v(1:2));
    t = cell_sub(T, I(1:2));
    % trimed tree
    T = cell_sub(T, I(3:length(I)) );
    p = v(3:length(v));
    % add a new node with the corresponding probability
    p(length(p)+1) = q;
    T = cell_set(T, length(p), t);
end

clf;
plot_hufftree(T);

C = huffman_gencode(T);

for i=1:size(C,1)
    disp(strcat(['Code of token ' num2str(i) ' = ' num2str( cell_get(C,i) )]));
end

n = 1024;

x = rand_discr(h, n);
x = x(:);

exo1()

%% Insert your code here.

e = - sum( h .* log2( max(h,1e-20) ) );
disp( strcat(['Entropy bound = ' num2str(n*e) '.']) );
disp( strcat(['Huffman code  = ' num2str(length(y)) '.']) );

t = cell_get(T,1);

x1 = [];

y1 = y;
while not(isempty(y1))
    % go down in the tree
    if y1(1)==0
        t = cell_get(t,1);
    else
        t = cell_get(t,2);
    end
    % remove the symbol from the stream buffer
    y1(1) = [];
    if not(iscell(t))
        % we are on a leaf of the tree: output symbol
        x1 = [x1 t];
        t = cell_get(T,1);
    end
end
x1 = x1(:);

err = norm(x-x1);
disp( strcat(['Error (should be 0)=' num2str(err) '.']) );

t = .12;

h = [t; 1-t];

n = 4096*2;
x = (rand(n,1)>t)+1;

q = 3;

m = 2;

n1 = ceil(n/q)*q;

x1 = x;
x1(length(x1)+1:n1) = 1;
x1 = reshape(x1,[q n1/q]);
[Y,X] = meshgrid(1:n1/q,0:q-1);
x1 = sum( (x1-1) .* (m.^X), 1 )' + 1;

H = h; 
for i=1:q-1
    Hold = H;
    H = [];
    for i=1:length(h)
        H = [H; Hold*h(i)];
    end
end

H = h;
for i=1:q-1
    H = kron(H,h);
end

exo2()

%% Insert your code here.

p = 0.1;

n = 512;

x = (rand(n,1)>p)+1;

h = [p 1-p];

y = perform_arith_fixed(x,h);

x1 = perform_arith_fixed(y,h,n);

disp(strcat(['Decoding error (should be 0)=' num2str(norm(x-x1)) '.']));

exo3()

%% Insert your code here.

n = 4096;

q = 10;
h = 1:q; h = h/sum(h);

x = rand_discr(h, n);

h1 = hist(x, 1:q)/n;
clf;
subplot(2,1,1); 
bar(h); axis('tight');
set_graphic_sizes([], 20);
title('True distribution');
subplot(2,1,2);
bar(h1); axis('tight');
set_graphic_sizes([], 20);
title('Empirical distribution');

warning off
exo4()
warning on

%% Insert your code here.
