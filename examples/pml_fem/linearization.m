%% Clear the workspace and the command window
clear;
clc;

%% Load the matrices
load('A0.mat', 'A0');
load('A1.mat', 'A1');
load('A2.mat', 'A2');

%% Linearization
% There are three linearizations (refer to the note)
n = size(A0, 1);
I = speye(n);

%% Linearization 1
% A = [A1, A0; -I, sparse(n, n)];
% B = [-A2, sparse(n, n); sparse(n, n), -I];

%% Linearization 2
A = [A1, -I; A0, sparse(n, n)];
B = [-A2, sparse(n, n); sparse(n, n), -I];

%% Liearization Kirsch & Zhang
% A = [A0, sparse(n, n); sparse(n, n), I];
% B = [-A1, -A2; I, sparse(n, n)];

%% Solve the generalized eigenvlaue problem
tic;
d = eigs(A, B, 200, 0);
toc;

%% Plot the eigenvalues

% Use LaTeX by default for all labels
set(groot,'defaultTextInterpreter','latex')
scatter(real(d), imag(d), 150, "red", ".");

% Set axis limit
axis([-0.5, 0.5, -0.5, 0.5]);
xlabel('$\Re(\alpha)$');
ylabel('$\Im(\alpha)$');
grid on;