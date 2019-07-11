function online_test()

clc, clear all, close all

training_output = 'training_output';
load_A = load(strcat(training_output,'/model_2_invariant_matrix_A.mat'));
load_B = load(strcat(training_output,'/model_2_invariant_matrix_B.mat'));

% load_matrices
A = load_A.A;
B = load_B.B;

end