clc, close all

svd_plots = 1;
sol_plots = 1;

%% LOAD FILES AND VARIABLES
load_sys_vars           = load(strcat(training_folder,'/model_1_sys_vars.mat'));
load_solution           = load(strcat(training_folder,'/model_1_solution_snapshots.mat'));
load_residual_snapshots = load('gqe_figs/residual_snapshots.mat');
residual_snapshots      = load_residual_snapshots.residual_snapshots; 
sys_vars                = load_sys_vars.sys_vars;
load_ul_ktr150 = load('gqe_figs/ul_ktr150.mat');
ul_ktr150 = load_ul_ktr150.ul_hdm;
ul_ktr150(:,3) = ul_ktr150(:,3) + 20;
load_ul_hdm = load('gqe_figs/ul_hdm.mat');
ul_hdm = load_ul_hdm.ul;
ul_hdm(:,3) = ul_hdm(:,3) + 20;
mesh_input               = sys_vars.mesh_input;


if svd_plots
    [Phi_s, Sigma_s, ~] = svd(load_solution.solution_snapshots, 'econ');
    [Phi_r, Sigma_r, ~] = svd(residual_snapshots, 'econ');

    podMode_s  = 1:size(Sigma_s,1);
    singVals_s = diag(Sigma_s);
    podMode_r  = 1:size(Sigma_r,1);
    singVals_r = diag(Sigma_r);

    % Find where to truncate given an error tolerance, if requested
    % sol
    error_s = zeros(size(Sigma_s,1),1);
    den_s = 0;
    tol=1e-6;
    for ii = 1 : size(Phi_s,2)
        den_s = den_s + Sigma_s(ii,ii)^2;
    end
    num_s = 0;
    errorSatisfied = false;
    for ii = 1 : size(Phi_s,2)
        num_s = num_s + Sigma_s(ii,ii)^2;
        error = 1 - num_s / den_s;
        error_s(ii) = error * 100;
        if error < tol && ~errorSatisfied
            ktrPossible = ii;
            errorSatisfied = true;
        end     
    end

    % res
    error_r = zeros(size(Sigma_r,1),1);
    den_r = 0;
    tol=1e-6;
    for ii = 1 : size(Phi_r,2)
        den_r = den_r + Sigma_r(ii,ii)^2;
    end
    num_r = 0;
    errorSatisfied = false;
    for ii = 1 : size(Phi_r,2)
        num_r = num_r + Sigma_r(ii,ii)^2;
        error = 1 - num_r / den_r;
        error_r(ii) = error * 100;
        if error < tol && ~errorSatisfied
            ktrPossible = ii;
            errorSatisfied = true;
        end     
    end

    %% SVD plots
    %% U plots
    figure(), hold on
    ymax = 1.05*singVals_s(1);
    semilogx(podMode_s,singVals_s,'k','LineWidth',2.5)
    plot(150*[1,1],[0,ymax],'r','LineWidth',1.5)
    xlabel('POD mode number','FontSize',25)
    ylabel('Singular value','FontSize',25)
    title('Singular values of snapshot solution matrix','FontSize',25)
    axis([0 1e3 0 ymax])
    set(gca, 'XScale', 'log');
    box on
    print('gqe_figs/therel_sing_vals_u.eps','-depsc')

    % Plot error as a function of number of basis modes
    figure(), hold on
    ymax = 100;
    semilogx(podMode_s,error_s,'k','LineWidth',2.5)
    plot(150*[1,1],[0,ymax],'r','LineWidth',1.5)
    xlabel('POD mode number','FontSize',25)
    ylabel('Relative percent error','FontSize',25)
    title('Relative percent error of solution ROM','FontSize',25)
    axis([0 1e3 0 ymax])
    set(gca, 'XScale', 'log');
    box on
    print('gqe_figs/therel_err_u.eps','-depsc')

    %% R plots
    figure(), hold on
    ymax = 1.05*singVals_r(1);
    semilogx(podMode_r,singVals_r,'k','LineWidth',2.5)
    plot(150*[1,1],[0,ymax],'r','LineWidth',1.5)
    xlabel('POD mode number','FontSize',25)
    ylabel('Singular value','FontSize',25)
    title('Singular values of snapshot residual matrix','FontSize',25)
    axis([0 1e3 0 ymax])
    set(gca, 'XScale', 'log');
    box on
    print('gqe_figs/therel_sing_vals_r.eps','-depsc')

    % Plot error as a function of number of basis modes
    figure(), hold on
    ymax = 100;
    semilogx(podMode_r,error_r,'k','LineWidth',2.5)
    plot(150*[1,1],[0,ymax],'r','LineWidth',1.5)
    xlabel('POD mode number','FontSize',25)
    ylabel('Relative percent error','FontSize',25)
    title('Relative percent error of residual ROM','FontSize',25)
    axis([0 1e3 0 ymax])
    set(gca, 'XScale', 'log');
    box on
    print('gqe_figs/therel_err_r.eps','-depsc')
end

if sol_plots
    %% MESH GEOMETRY
    [sys_vars,xl_hdm,connect_list,element_dof,dof_list] = mesh(sys_vars,mesh_input);
    nXele = 10;
    nYele = 10;
    xDim = 0.1;
    yDim = 0.1;

    diagonalIndices = zeros(nYele,1);
    diagonalIndices(1) = 1;
    for ii = 2 : nYele+1
        diagonalIndices(ii) = diagonalIndices(ii-1) + nXele + 2;
    end
    diagonalDistance = (0:nYele)/nYele * sqrt(xDim^2 + yDim^2);% - sqrt(xDim^2 + yDim^2)/nYele;

    % with deim. best line
    figure(), hold on
    plot(diagonalDistance,ul_ktr150(diagonalIndices,3),'k','LineWidth',2)
    plot(diagonalDistance,full(ul_hdm(diagonalIndices,3)),'r--','LineWidth',2)

    axis([0 sqrt(2)*0.1 20 28])
    xlabel('Distance along diagonal','FontSize',25)
    ylabel('Temperature (C)','FontSize',25)
    lgd = legend('HDM','ns=100','Location','northwest');
    lgd.FontSize = 15;
    legend('boxoff')
    theTitle = strcat('Temperature along diagonal');
    title(theTitle,'FontSize',25)
    box on
    print('gqe_figs/temp_along_diag.eps','-depsc')
    
    error = norm(full(ul_hdm)-ul_ktr150) / norm(full(ul_hdm)) * 100


end
