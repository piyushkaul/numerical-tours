
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('toolbox_graph')
addpath('solutions/shapes_4_shape_matching')

n = 250;
name = {'square' 'disk'};
name = {'centaur1' 'centaur2'};
M = rescale( load_image(name,n) );

M = perform_blurring(M,3);
M = {double(M{1}>.5) double(M{2}>.5)};

clf;
imageplot(M, {'Shape 1', 'Shape 2'});

for i=1:2
    bound0{i} = compute_shape_boundary(M{i});
end

nbound = 500;
bound = {};
for i=1:2
    nbound0 = size(bound0{i},2);
    t = linspace(1,nbound0,nbound+1); t(end) = [];
    bound{i}(1,:) = interp1( bound0{i}(1,:), t );
    bound{i}(2,:) = interp1( bound0{i}(2,:), t );
end

for i=1:2
    [tmp,pt] = min(bound{i}(1,:));
    bound{i} = [bound{i}(:,pt:end) bound{i}(:,1:pt-1)];
end

clf;
for i=1:2
    b = bound{i};
    subplot(1,2,i);
    hold on;
    imageplot(M{i});
    hh = plot(b(2,:), b(1,:)); 
    set(hh, 'LineWidth', 2);
    hh = plot(b(2,1), b(1,1), '.r'); 
    set(hh, 'MarkerSize', 20);
    axis('ij');
end

nrad = 10;
rlist = linspace(3,20,nrad);

k = 4;
r = rlist(k);

x = -ceil(r):ceil(r);
[b,a] = meshgrid(x,x);
h = double( a.^2+b.^2<=r^2 );
h = h/sum(h(:));

Mh = perform_convolution(M,h);

clf;
imageplot(Mh,{'Blurred 1', 'Blurred 2'});

[Y,X] = meshgrid(1:n,1:n);
D = {};
for i=1:2
    D{i}(k,:) = interp2(Y,X,Mh{i},bound{i}(2,:), bound{i}(1,:));
end

exo1()

%% Insert your code here.

sel = ceil( [1 90 400] * nbound/500 );
col = {'r', 'g', 'b', 'k'};
clf; hold on;
imageplot(M{1});
for i=1:length(sel)
    h = plot(bound{1}(2,sel(i)), bound{1}(1,sel(i)), [col{i} '.']);
    set(h, 'MarkerSize', 40);
end
axis('ij');

clf;
col = {'r', 'g', 'b', 'k'};
a = mean(mean(D{1}(:,sel)));
for i=1:length(sel)
    subplot(length(sel),1,i);
    h = bar(D{1}(:,sel(i))-a, col{i}); axis('tight');
    set(gca, 'XTickLabel', []);
end

exo2()

%% Insert your code here.

exo3()

%% Insert your code here.

hold on;
imageplot(C);
h = plot(gpath(2,:), gpath(1,:), 'r');
set(h, 'LineWidth', 2);
axis('xy');

clf;
col = {'r' 'b'};
clf; hold on;
delta = [0 n/2];
for i=1:2
    h = plot3(bound{i}(1,[1:end 1]), bound{i}(2,[1:end 1]), delta(i)*ones(nbound+1,1), col{i});
    set(h, 'LineWidth', 2);
end
t = round(linspace(length(gpath)/2.1,length(gpath),60)); 
for i=1:length(t)
    i1 = round(gpath(1,t(i)));
    i2 = round(gpath(2,t(i)));
    a = bound{1}(:,i1);
    b = bound{2}(:,i2);
    h = plot3( [a(1) b(1)], [a(2) b(2)], delta, 'k-'  );
    set(h, 'LineWidth', 1);
end
view(-10,40);
axis('equal'); axis('off');
zoom(1.3);
