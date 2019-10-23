
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/introduction_6_elementary_fr')

n = 256;
name = 'hibiscus';
f = load_image(name, n);
f = rescale(sum(f,3));
clf;
imageplot(f);

selx = 19:24;
sely = 62:67;
clf;
image(f(selx,sely)*255); axis image; axis off;
disp(floor(255*f(selx,sely)));

sub = @(f,k)f(1:k:end,1:k:end);

clf;
imageplot(sub(f,4));

klist = [2 4 8 16];
clf;
for i=1:length(klist)
    k = klist(i);
    imageplot(clamp(sub(f,k)), ['1 ligne/colonne sur ' num2str(k)], 2, 2, i);
end

quant = @(f,q)(round(q*rescale(f,1e-3,1-1e-3)-1/2)+1/2)/q;
clf;
imageplot(quant(f,4), '4 niveaux de gris');

qlist = [16, 4, 3, 2];
clf;
for i=1:length(qlist)
    q = qlist(i);
    f1 = quant(f,q); f1(1)=0; f1(2)=1;
    imageplot(f1, [num2str(q) ' niveaux de gris' ], 2, 2, i);
end

name = 'boat';
f = rescale(load_image(name, n));
sigma = .08;
f = f + randn(n)*sigma;
clf;
imageplot(f);

filt_moy = @(f,k)perform_convolution(f,ones(2*k+1)/(2*k+1)^2,'sym');
clf;
imageplot(clamp(f), 'Image bruit e', 1, 2, 1);
imageplot(clamp(filt_moy(f,1)), 'Image moyenn e', 1, 2, 2);

klist = [1 2 3 4];
clf;
for i=1:length(klist)
    k = klist(i);
    f1 = filt_moy(f,k);
    imageplot(clamp(f1), ['Moyenne de ' num2str((2*k+1)^2) ' pixels'], 2, 2, i);
end

filt_med = @(f,k)perform_median_filtering(f,k);
clf;
imageplot(clamp(filt_moy(f,1)), 'Moyenne de 9 nombres', 1, 2, 1);
imageplot(clamp(filt_med(f,1)), 'M diane de 9 nombres', 1, 2, 2);

klist = [1 2 3 4];
clf;
for i=1:length(klist)
    k = klist(i);
    f1 = filt_med(f,k);
    imageplot(clamp(f1), ['M diane de ' num2str((2*k+1)^2) ' pixels'], 2, 2, i);
end

n = 256;
name = 'hibiscus';
f = load_image(name, n);
f = rescale(sum(f,3));

s1 = [1 1:n-1]; s2 = [2:n n];
edge = @(f)sqrt( ( f(s1,:) - f(s2,:) ).^2 + ( f(:,s1) - f(:,s2) ).^2 );

clf;
imageplot(f, 'Image', 1,2,1);
imageplot(edge(f), 'Carte de l', 1,2,2);

name = 'hibiscus';
f = rescale( load_image(name,n) );

f1 = cat(3, f(:,:,1), zeros(n), zeros(n));
f2 = cat(3, zeros(n), f(:,:,2), zeros(n));
f3 = cat(3, zeros(n), zeros(n), f(:,:,3));

clf;
imageplot({f f1 f2 f3}, ...
        { 'Image couleur' 'Canal rouge' 'Canal vert' 'Canal bleu'}, 2, 2);

clf;
imageplot({f sum(f,3)}, {'Couleur' 'Luminance'});

g = 1-f;
f1 = cat(3, f(:,:,1),     f(:,:,2)*0+1, f(:,:,3)*0+1);
f2 = cat(3, f(:,:,1)*0+1, f(:,:,2)    , f(:,:,3)*0+1);
f3 = cat(3, f(:,:,1)*0+1, f(:,:,2)*0+1, f(:,:,3));

clf;
imageplot({f f1 f2 f3}, ...
        { 'Image couleur' 'Canal cyan' 'Canal magenta' 'Canal jaune'}, 2, 2);

name = 'hibiscus';
f = rescale( load_image(name,n) );
f = rescale(sum(f,3));

clf;
imageplot(-f);

clf;
imageplot(f.^2, 'Carr ');

clf;
imageplot(sqrt(f), 'Remplacement de a par sqrt(a)');

name = 'hibiscus';
f = rescale( load_image(name,n) );

m = @(f)repmat(mean(f,3), [1 1 3]);
contrast = @(f,gamma)clamp(m(f).^gamma + f-m(f));

gamma_list = [.5 .75 1 1.5 2 3];
clf;
for i=1:length(gamma_list)
    subplot(2,3,i);
    image(contrast(f,gamma_list(i))); axis image; axis off;
    title(['\gamma=' num2str(gamma_list(i))]);
    colormap jet(256);
end

A = rescale( load_image('flowers',512) );
B = permute(A, [2 1 3]);
clf;
imageplot({A B}, {'Image A' 'Image B'}, 1,2,1);

C = A;
C = C(end:-1:1,:,:); C = permute(C, [2 1 3]); 
clf;
imageplot({A C}, {'Image A' 'Image C'}, 1,2,1);

A = rescale( load_image('flowers',512) );
B = rescale( load_image('hibiscus',512) );
clf;
imageplot({A B}, {'Image A', 'Image B'});

p = 6;
t = linspace(0,1,p);
clf;
for i=1:p
    imageplot(t(i)*A+(1-t(i))*B, ['t=' num2str(t(i), 2)], 2,p/2,i);
end
