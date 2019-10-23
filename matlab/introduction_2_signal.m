
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/introduction_2_signal')

n = 512;

f = load_signal('Piece-Regular', n); % signal of size n

f = f(:);

f = rescale(f);

clf;
plot(1:n, f);
axis('tight');
title('My title'); % title
set_label('variable x', 'variable y'); % axis

subplot(2, 2, 1);
plot(f); axis('tight');

subplot(2, 2, 4);
plot(f.^2); axis('tight');

clf;
plot(1:n, [f f.^2]');
legend('signal', 'signal^2');
