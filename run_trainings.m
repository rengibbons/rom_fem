%% CLEAR WORKSPACE
clc, close all
format short

%% CHOOSE WHICH MODEL(S) TO RUN
run_mod_1   = 0;
run_mod_2   = 0;
run_mod_3   = 1;

%% FLAGS
flags.debug                      = 0;
flags.enable_numerical_tangent   = 0;
flags.write_tec                  = 0;
flags.plot_load                  = 0;
flags.plot_initial_configuration = 0;
flags.plot_each_load_step        = 0;
flags.plot_final_configuration   = 1;
flags.save_solution_vectors      = 1;

%% INPUT FILE

% ktr = 27; % used for gqe neohookean elastic test
% input_file      = 'input_file_elastic_test_neoh_10x10';
% training_folder = 'training_folder_10x10';

% ktr = 20; % good for ktr = 20
% input_file      = 'input_file_thermo_elastic_test_04x04_ther_loading';
% training_folder = 'training_folder_thermo_elastic_04x04_ther_loading';

% ktr = 35; % goodish for ktr = 35
% input_file      = 'input_file_thermo_elastic_test_04x04_ther_mech_loading';
% training_folder = 'training_folder_thermo_elastic_04x04_ther_mech_loading';

% ktr = 50; % good for ktr = 50
% input_file      = 'input_file_thermo_elastic_test_10x10_ther_mech_loading';
% training_folder = 'training_folder_thermo_elastic_10x10_ther_mech_loading';

ktr = 200; % good for ktr = 150
input_file      = 'input_file_thermo_elastic_10x10_moving_heat_src';
% training_folder = 'train_ther_el_10x10_moving_heat_src_20x20_samples_alpha_1e-3';
% training_folder = 'train_ther_el_10x10_moving_heat_src_5x5_samples_alpha_1e-3';
% training_folder = 'train_ther_el_10x10_moving_heat_src_10x10_samples_alpha_1e-4';
training_folder = 'train_ther_el_10x10_moving_heat_src_5x5_samples_alpha_1e-7';
% training_folder = strcmp(training_folder,'_ktr_',num2str(ktr))

% ktr = 300; % good for ktr = ???
% input_file      = 'input_file_thermo_elastic_16x16_moving_heat_src';
% training_folder = 'train_ther_el_16x16_moving_heat_src_5x5_samples_alpha_1e-7';

% ktr = 500; % good for ktr = ???
% input_file      = 'input_file_thermo_elastic_20x20_moving_heat_src';
% training_folder = 'train_ther_el_20x20_moving_heat_src_10x10_samples_alpha_1e-4';

% ktr = 250; % good for ktr = 
% input_file      = 'input_file_thermo_elastic_64x64_moving_heat_src';
% training_folder = 'train_ther_el_64x64_moving_heat_src_5x5_samples_alpha_1e-4';

% ktr = 250; % good for ktr = 
% input_file      = 'input_file_thermo_elastic_32x32_moving_heat_src';
% training_folder = 'train_ther_el_32x32_moving_heat_src_5x5_samples_alpha_1e-4';

%%
training_params.nu_training    = [0 0.1 0.2 0.3 0.4 0.45]';
training_params.alpha_training = [1e-3 5e-3 1e-2]'; % not used. never tested

x0 = [0.09 0.09];
xf = [0.01 0.01];
n_dim = size(x0,2);

n_disc = 5;  % Model I: 10 mins
% n_disc = 10;  % Model I: ???
% n_disc = 20; % Model I: 120 mins
x_disc = sort(linspace(x0(1),xf(1),n_disc));

[x,y] = meshgrid(x_disc);
training_params.heat_src_training = zeros(size(x,1)*size(x,2),n_dim);
training_params.heat_src_training(:,1) = x(:);
training_params.heat_src_training(:,2) = y(:);

%% RUN MODELS
if run_mod_1
    run_model_1(input_file,training_folder,training_params,flags)
end

% flags.plot_final_configuration   = 0;

if run_mod_2
    run_model_2(training_folder,ktr,training_params,flags)
end

flags.plot_each_load_step        = 0;
flags.plot_final_configuration   = 1;

if run_mod_3
    run_model_3(training_folder,flags)
end

%%
% ktr = 2;
% input_file      = 'patch_test';
% training_folder = 'training_folder_patch_test';

% ktr = 15;
% input_file      = 'input_file_elastic_test_neoh_03x03'; % good for ktr = 15
% training_folder = 'training_folder_03x03';

% ktr = 20;
% input_file      = 'input_file_elastic_test_neoh_04x04'; % good for ktr = 20
% training_folder = 'training_folder_04x04';

% ktr = 30;
% input_file      = 'input_file_elastic_test_neoh_05x05'; % good for ktr = 30
% training_folder = 'training_folder_05x05';

% ktr = 50;
% input_file      = 'input_file_elastic_test_neoh_07x07';
% training_folder = 'training_folder_07x07';

