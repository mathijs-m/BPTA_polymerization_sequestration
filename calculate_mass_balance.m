function [c_BPTA, c_additive, COUNT] = calculate_mass_balance(K_pol_BPTA, sig_BPTA, K_pol_additive, sig_additive, K_seq, c_BPTA, c_additive, sequestrating_stoichiometry)
% Set initial errors equal to function values at origin.
err_c_BPTA = calculate_mass_balance_BPTA(0, c_BPTA,0, K_pol_BPTA, sig_BPTA, K_seq, sequestrating_stoichiometry);
err_c_additive = calculate_mass_balance_additive(0, c_additive,0, K_pol_additive, sig_additive, K_seq, sequestrating_stoichiometry);

% Set an initial starting point
c_BPTA_test = 0;
c_additive_test = 0;

% Set the reference point very far. This is just to enter the while loop
c_BPTA_test_old = 1e99;
c_additive_test_old = 1e99;

tol = 1e-8;

err_c_BPTA = [];
err_c_additive = [];
optimized_c_BPTA = [];
optimized_c_additive = [];
COUNT = 0;

% Perform the nested binary search to find values for c_BPTA and c_additive
% for which both mass-balance equations are satisfied. Iterate until the
% tolerances are met.
while ~all([abs(err_c_BPTA)<0.1*tol,...
            abs(err_c_additive)<tol,...
            abs(c_BPTA_test-c_BPTA_test_old)<1e-5*c_BPTA,...
            abs(c_additive_test-c_additive_test_old)<c_additive*1e-5])    
    
    % Define the search boundaries
    c_BPTA_guess = [0, min([1/K_pol_BPTA,c_BPTA])-1e-30];
    
    % Save the solution from the previous iteration
    c_BPTA_test_old = c_BPTA_test;
    c_additive_test_old = c_additive_test;%
    
    c_BPTA_test = Asolver(c_BPTA_guess,c_BPTA,c_additive_test, K_pol_BPTA, sig_BPTA, K_seq, model,tol);
    optimized_c_BPTA = [optimized_c_BPTA,c_BPTA_test];
    
    Bguess = [0, min([1/K_pol_additive,c_additive])-1e-20];
    c_additive_test = Bsolver(Bguess,c_additive,c_BPTA_test, K_pol_additive, sig_additive, K_seq, model, tol);
    optimized_c_additive = [optimized_c_additive, c_additive_test];
    
    % Calculate the errors at the current [Mtest,Stest] coordinates
    err_c_BPTA = Acalculator(c_BPTA_test,c_BPTA,c_additive_test,K_pol_BPTA, sig_BPTA, K_seq, model);
    err_c_additive = Bcalculator(c_additive_test,c_additive,c_BPTA_test,K_pol_additive, sig_additive, K_seq, model);

    COUNT = COUNT + 1; % Master counter
    err_c_BPTA = [err_c_BPTA, err_c_BPTA];
    err_c_additive = [err_c_additive,err_c_additive];   

    if COUNT > 100000 && any([c_BPTA,c_additive] == 0)
       % Break out of the loop after a certain limit or the concentration
       % are found to be 0, i.e., the system is aggregating too strongly.
       display('Count overflow.')       
       break
    end
end
c_BPTA = c_BPTA_test;
c_additive = c_additive_test;
end

function Aopt = Asolver(Aguess,c_a,Btest, Ka, sig_a, Kseq, model,tol)
% Calculate the error in M
errA_L = Acalculator(Aguess(1),c_a,Btest, Ka, sig_a, Kseq, model);
errA_R = Acalculator(Aguess(2),c_a,Btest, Ka, sig_a, Kseq, model);

count = 0;
while any([errA_L>0.0001*tol, ((Aguess(2)-Aguess(1))>1e-15)])
    xA = (Aguess(1) + Aguess(2))/2;
    errA = Acalculator(xA,c_a,Btest, Ka, sig_a, Kseq, model);

    if errA > 0
        Aguess(1) = xA;
        errA_L = errA;
    else
        Aguess(2) = xA;
        errA_R = errA;
    end
    count = count + 1;
    if count > 1e3
        display('Counter overflow for A')
        break
    end
end

Aopt = Aguess(1);
end

function Bopt = Bsolver(Bguess,c_b,Atest, Kb, sig_b, Kseq, model, tol)
    count = 0; % Reset the counter

    % Calculate the error in S
    errB_L = Bcalculator(Bguess(1),c_b,Atest, Kb, sig_b, Kseq, model);
    errB_R = Bcalculator(Bguess(2),c_b,Atest, Kb, sig_b, Kseq, model);
    
    % Do the binary search in S
    while any([errB_L>0.0001*tol, ((Bguess(2)-Bguess(1))>1e-10)])    
        xB = (Bguess(1) + Bguess(2))/2;
        errB = Bcalculator(xB,c_b,Atest,Kb, sig_b, Kseq, model);
        
        if errB > 0
            Bguess(1) = xB;
            errB_L = errB;
        else
            Bguess(2) = xB;
            errB_R = errB;
        end
        
        count = count + 1;
        if count > 1e2
            error('Counter overflow for B')
            break
        end
    end
    % Save the solution
    Bopt = Bguess(1);
end



%% 
%ACALCULATOR%
function err_BPTA = calculate_mass_balance_BPTA(c_BPTA,c_BPTA_tot,c_additive, K_pol_BPTA, sig_BPTA, K_seq, sequestrating_stoichiometry)
BPTA_pol = (1-sig_BPTA).*c_BPTA+sig_BPTA.*c_BPTA./(1-K_pol_BPTA.*c_BPTA).^2;
if sequestrating_stoichiometry == 1
    BPTA_seq = K_seq.*c_BPTA.*c_additive;
elseif sequestrating_stoichiometry == 2
    BPTA_seq = K_seq.*c_BPTA.*c_additive+2*K_seq.^2*c_BPTA.^2*c_additive;
elseif sequestrating_stoichiometry == 3
    BPTA_seq = K_seq.*c_BPTA.*c_additive+2*K_seq.^2*c_BPTA.^2*c_additive + 3.*K_seq.^3.*c_BPTA.^3.*c_additive;
end
err_BPTA = c_BPTA_tot - BPTA_pol - BPTA_seq;

end

%BCALCULATOR%
function y = calculate_mass_balance_additive(c_additive,c_additive_tot,c_BPTA,K_pol_additive, sig_additive, K_seq, sequestrating_stoichiometry)
additive_pol = (1-sig_additive).*c_additive+sig_additive.*c_additive./(1-K_pol_additive.*c_additive).^2;
if sequestrating_stoichiometry == 1
    additive_seq = K_seq.*c_BPTA.*c_additive;
elseif sequestrating_stoichiometry == 2
    additive_seq = K_seq.*c_BPTA.*c_additive+K_seq.^2*c_BPTA.^2*c_additive;
elseif sequestrating_stoichiometry == 3
    additive_seq = K_seq.*c_BPTA.*c_additive+2*K_seq.^2*c_BPTA.^2*c_additive + K_seq.^3.*c_BPTA.^3.*c_additive;
end
y = c_additive_tot - additive_pol - additive_seq;

end

