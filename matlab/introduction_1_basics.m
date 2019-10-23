
addpath('toolbox_signal')
addpath('toolbox_general')
addpath('solutions/introduction_1_basics')

a = 1; a = 2+1i; % real and complex numbers
b = [1 2 3 4]; % row vector
c = [1; 2; 3; 4]; % column vector
d = 1:2:7; % here one has d=[1 3 5 7]

size(d)

d(1)

d(1:2)

A = eye(4,4); 
B = ones(4,4);
C = rand(4,4);

c = b';

D = C*A
D = C.*A

E = A./C; % division
E = sin(A); % sinus is applied to each entry
E = abs(A + 1i*C); % modulus of each entry

b = sort(b); % sort values
b = b .* (b>2); % set to zeros (threshold) the values below 2
b(3) = []; % suppress the 3rd entry of a vector
B = [b; b]; % create a matrix of size 2x4
c = B(:,2); % to access 2nd column

b(end-2:end) = 1; % to access the last entries

b = b(end:-1:1); % reverse a vector

disp('Hello'); % display a text
x = 1.23456;
disp( sprintf('Value of x=%.2f', x) ); % print a values with 2 digits
A(A==Inf) = 3; % replace Inf values by 3
A(:); % flatten a matrix into a column vector
max(A(:)); % max of a matrix

C = C .* (abs(C)>.3);

for i=1:3 % repeat the loop for i=1, i=2, i=3
    disp(i); % make here something
end
i = 3;
while i>0 % while syntax
    disp(i); % do smth
    i = i-1;
end

n = 256; % size of the image
M = load_image('lena', n);
clf;
imageplot(M);

clf;
imageplot(M(1:50,1:50), 'Zoom', 1,2,1);
imageplot(-M, 'Reversed contrast', 1,2,2);
