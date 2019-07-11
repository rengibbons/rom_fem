%% CLEAR WORKSPACE
clc, close all

%% INPUT FILE
% input_file = 'input_file_elastic_test_neoh_04x04';
% training_folder = 'training_folder_04x04';

% input_file_elastic_test_neoh_10x10
% training_folder = 'training_folder_10x10';

input_file    = 'input_file_thermo_elastic_test_04x04';
training_folder = 'training_folder_thermo_elastic_04x04';

% input_file_elastic_test_neoh
% input_file_thermo_elastic_test
% input_file_moving_heat_source

run_model_1(input_file,training_folder)
run_model_2(training_folder)