%%This program fits competition between a cooperative and a isodesmic pathway
%to a series of melting curves at different concentration in a single
%solvent, using 1 CD and 1 UV channel.

clear all
close all
clc

%% load and extract data
%select the folder containing this script. Ensure the experimental data is
%saved in a folder named 'data', Similarly, create two folders named 'figures'
%and 'results' respectively. These folders will hold the files outputted by the
%system. 

%extract filenames from 'data'-folder
data = [];
conc = [];
datasets = {'DataPhePyTA.txt','DataPhePyTA_40.txt','SelectedDataPhePyTA.txt'}
for dataset_index = 1:3
	dataset = datasets{dataset_index};
	data = dlmread(dataset);
	tic
	for sequestrating_stoichiometry = 1:3
        mkdir([dataset,'_',num2str(sequestrating_stoichiometry)]);
        %% input parameters and constants
        R = 8.3145; %gas constant J/mole/K
        PL = 1;     %cm
        total_samples = 200; %number of parameter sets

        E_mono = 0; % Molar ellipticity monomers
        E_seq = 0; % Molar ellipticity sequestrated
        E_pol_BPTA = 63.2586/30e-6; % Molar ellipticity polymer (Taken from data at 26.67 ?C)
        E_pol_additive = 40/300e-6; % Molar ellipticity polymer (Taken from data at 26.67 ?C)

        H_BPTA = -86e3;
        S_BPTA = -155;
        NP_BPTA = 15.4e3;
        H_additive = -94e3;  
        S_additive = -220;
        NP_additive = NP_BPTA;

        c_tot = 30e-6;

        %% generate fit parameters
        %input sample range
        H_max = -70e3; %estimated Gibbs free energy at 0%
        H_min = -100e3;

        S_min = -100; % m value for the aggregation at low vol%
        S_max = -200; 

        %Generate sample space us latin hypercube sampling which generates evenly
        %space values between 0 and 1
        intitial_parameters = lhsdesign(total_samples,2);

        intitial_parameters(:,1) = H_min + (H_max-H_min).*intitial_parameters(:,1);
        intitial_parameters(:,2) = S_min + (S_max-S_min).*intitial_parameters(:,2);

        %generate enthalpies from estimated Gibbs free energies and entropies
        axis_labels = {'H', 'S'};

        %define fitting settings
        constants = [R, PL, E_mono, E_seq, E_pol_BPTA, E_pol_additive, c_tot, H_BPTA, S_BPTA, NP_BPTA, H_additive, S_additive, NP_additive, sequestrating_stoichiometry];
        max_iter = 50;
        best_resnorm = 10e99;
        check = [];
        all_resnorms = [];
        all_optimized_parameters= [];
        all_residuals = [];
        all_jacobians = cell(1,1);
        Freport = [];
        Sreport = [];
        lowerbound = [H_min*ones(total_samples,1), S_min.*ones(total_samples,1)];
        upperbound = [H_max*ones(total_samples,1), S_max.*ones(total_samples,1)];
        curent_sample = 1;
        j = 1;

        %% Fitting
        while curent_sample < total_samples
            test_parameters = intitial_parameters(curent_sample,:);%load current parameters        

            %define optimisation options
            output_fun= @(x,optimValues,state,varargin)interuptFun(x,optimValues,state,varargin,best_resnorm,max_iter);
            options = optimset('MaxIter', 200, 'Display', 'iter', 'MaxFunEvals', 2000, ...
            'TolFun', 1e-15, 'TolX', 1e-15, 'Algorithm', 'levenberg-marquardt', 'OutputFcn', output_fun, 'UseParallel', true);

            lb = [];%lowerbound(n,:);
            ub = [];%upperbound(n,:);

            try
                [par_fin,resnorm,residual,exitflag,output,lambda,Jacobian] =... 
                lsqnonlin(@calculate_cost_function,test_parameters, lb,ub, options, ...
                data, constants);

                %save characteristics of current best fit
                if resnorm < best_resnorm
                    best_resnorm = resnorm;
                    best_index = curent_sample;
                    check = [check, curent_sample];
                    best_jacobian = Jacobian;
                    best_residual = residual;
                    best_exit_flag = exitflag;
                end
                clc
                display(strcat('Progress: ', num2str(curent_sample/total_samples*100), '%'));
                display(['Resnorm: ',num2str(best_resnorm)])
                toc

                %save all succesfully optimised parameter sets
                all_resnorms = [all_resnorms; resnorm];
                Sreport = [Sreport; test_parameters];
                if ~isrow(par_fin)
                    all_optimized_parameters = [all_optimized_parameters; par_fin'];
                else
                    all_optimized_parameters = [all_optimized_parameters; par_fin];
                end
                all_residuals = [all_residuals; residual'];
                all_jacobians{1,curent_sample} = Jacobian; 

                curent_sample = curent_sample + 1;
                j = 1;

            catch
                Freport = [Freport; test_parameters];
                new_par = lhsdesign(1,2);

                newpar(:,1) = H_min + (H_max-H_min).*newpar(:,1);
                newpar(:,2) = S_min + (S_max-S_min).*newpar(:,2);

                intitial_parameters(curent_sample,:) = newpar;

                j = j + 1;
                if j > 1
                    j = 1;
                    curent_sample = curent_sample + 1;
                end
            end
        end

        %% Analysis of fit
        try
            best_finalized_parameters = all_optimized_parameters(best_index,:);
            best_residual = all_residuals(best_index,:);
        catch
            error('No fit was found')
        end
        storage.best_finalized_parameters = best_finalized_parameters;
        best_fit_indices = find(all_resnorms == best_resnorm);

        %Check for multiple distinct minima
        if length(best_fit_indices) > 1
            all_best_fits = all_optimized_parameters(best_fit_indices,:);
            differences_between_best_fits = all_best_fits - ones(size(all_best_fits),1)*best_finalized_parameters;
            if sum(sum(differences_between_best_fits)) ~= 0
                warning('Multiple minima detected for' + num2str(conc) + 'C');
            end
        end

        %Find all good fits
        good_indices = find(all_resnorms <= 1.01*best_resnorm);
        good_fits = all_optimized_parameters(good_indices,:);
        log_good_fits = sign(good_fits).*log10(abs(good_fits));
        storage.GoodFits = good_fits; 

        %Generate some figures
        if length(good_indices) == 1
            warning(['No boxplot is generated because there are no fits within 5% of the best fit']);
        elseif ~isempty(log_good_fits) && size(log_good_fits,1)>2
            figure;                       %Create a Boxplot of goodfits
            boxplot(log_good_fits, 'Labels', {axis_labels});
            title('Boxplot of logarithms of fit parameters within 1% of the best fit');
            savefig([dataset,'_',num2str(sequestrating_stoichiometry),'\2Pathways_box.fig'])
            close(gcf);
            figure;    %Create matrix plot which compares the optimized parameters.
            [S,AX,BigAx,H,HAx] = plotmatrix(good_fits);
            title(BigAx,'comparisons of fit parameters')
            for j = 1:size(good_fits,2)
                AX(j,1).YLabel.String = axis_labels(j);
                AX(size(good_fits,2),j).XLabel.String = axis_labels(j);
            end
            savefig([dataset,'_',num2str(sequestrating_stoichiometry),'\2Pathways_cor.fig']);
            close(gcf);
        else
            disp(['Only one distinct fit found for ', num2str(conc),', no boxplot generated']);
        end

        %Extract best fit parameters
        H_fit = best_finalized_parameters(1);
        S_fit = best_finalized_parameters(2);

        %Calculate bestfit curves
        uniquefracs = unique(data(:,1));
        for i = 1:length(uniquefracs)
            frac = uniquefracs(i);
            sim_T = linspace(0,100,100);
            simulation_data = [ones(1,100)*frac;sim_T]';
            [calculated_CD, calc_BPTA, calc_additive, calc_polymerized_BPTA, calc_polymerized_additive, calc_sequestrated_BPTA, calc_sequestrated_additive] = simulate_polymerization_sequestration(best_finalized_parameters, constants,simulation_data);
            figure(1);
            hold on
            plot(sim_T, calculated_CD, '-');
            xlabel('Temperature (C)'); ylabel('Absorbance (a.u.)');
            figure(2);
            hold on
            plot(sim_T, calc_BPTA,'r',sim_T, calc_additive,'m',sim_T,calc_polymerized_BPTA,'g',sim_T,calc_polymerized_additive,'y',sim_T,calc_sequestrated_BPTA,'b',sim_T,calc_sequestrated_BPTA,sim_T,(calc_BPTA+calc_polymerized_BPTA+calc_sequestrated_BPTA),'k',sim_T,calc_additive+calc_polymerized_additive+calc_sequestrated_additive);
            legend('MonomerA','MonomerB', 'PolymerA', 'PolymerB','SequestratedA', 'SequestratedB','TotalA','TotalB')
            xlabel('Temperature (?C)'); ylabel('Concentration');

            %save the simulated signal
            FileID = fopen(strcat([dataset,'_',num2str(sequestrating_stoichiometry),'\2Pathways_signal_c=_', num2str(frac), '.txt']),'wt');
            txt = sprintf('Fraction \t Absorbance signal');
            fprintf(FileID, '%s\t\r\n' ,txt);
            for j = 1:length(sim_T)
                DATA = [sim_T(j),calculated_CD(j)];
                fprintf(FileID, '%5f\t%5f\n', DATA);
            end

            %species
            FileID = fopen(strcat([dataset,'_',num2str(sequestrating_stoichiometry),'\PhePyTA_sequestration_ratio=', num2str(frac), '.txt']),'wt');
            txt = sprintf('Temperature \t MonomerA \t MonomerB \t PolymerA \t PolymerB \t SequestratedA \t SequestratedB \t TotalA \t TotalB');
            fprintf(FileID, '%s\t\r\n', txt);
            for j = 1:length(sim_T)
                DATA = [sim_T(j), calc_BPTA(j),calc_additive(j), calc_polymerized_BPTA(j), calc_polymerized_additive(j), calc_sequestrated_BPTA(j),calc_sequestrated_additive(j), (calc_BPTA(j) + calc_polymerized_BPTA(j) + calc_sequestrated_BPTA(j)),(calc_additive(j) + calc_polymerized_additive(j) + calc_sequestrated_additive(j))];
                fprintf(FileID, '%5e\t%5e\t%5e\t%5e\t%5e%5e\t%5e\t%5e\t\r\n', DATA);
            end

        end
    figure(1)
    plot(data(:,2),data(:,3),'o','MarkerSize',2)
    savefig(strcat([dataset,'_',num2str(sequestrating_stoichiometry),'\2Pathways_Fit_c=_', num2str(frac), '.fig']));
    figure(2)
    savefig(strcat([dataset,'_',num2str(sequestrating_stoichiometry),'\2Pathways_Speciation_c=_', num2str(frac), '.fig']));

    %perform statistics on the fit parameters
    df = -length(best_finalized_parameters); %number of degrees of freedom
    for i = 1:size(data,2)
        df = df + length(data(:,1));
    end
    alpha = 0.32; %one standard deviation
    try
        ci = nlparci(best_finalized_parameters, best_residual, best_jacobian, alpha); %calculates the 2*alpha% confidence interval, returns upper and lower boundaries per parameter
        tsd = tinv(1-alpha./2,df); 
        x_sd = (ci(:,2)-ci(:,1))./(2*tsd); % calculates the standard deviation for each parameter
    catch
        x_sd = zeros(size(best_finalized_parameters));
    end
    storage.x_sd = x_sd;
    storage.Jacob = best_jacobian; 
    storage.Constants = constants;

    save('storage_iso.mat', '-struct','storage');

    %export parameters
    FileID = fopen(strcat([dataset,'_',num2str(sequestrating_stoichiometry),'\2Pathways_parameters.txt']), 'wt');
    txt = sprintf('H= %5e\t+/-\t%5e\r\nS= %5e\t+/-\t%5e\r\n\nresnorm= %5e', ...
        H_fit, x_sd(1),S_fit,x_sd(2),best_resnorm);
    fprintf(FileID,'%s\t\r\n', txt);
    fclose(FileID);
    close all
    end
end