function [CDcalc,free_BPTA,free_additive,polymerized_BPTA, polymerized_additive,sequestrated_BPTA,sequestrated_additive] = simulate_polymerization_sequestration(par,constants,data)

H_seq = par(1);
S_seq = par(2);

R = constants(1);
PL = constants(2);

Emono = constants(3);
Eseq = constants(4);
Epol_a = constants(5);
Epol_b = constants(6);

c_tot = constants(7);
H_BPTA = constants(8);
S_BPTA = constants(9);
NP_BPTA = constants(10);
H_additive = constants(11);
S_additive = constants(12);
NP_additive = constants(13);
sequestrating_stoichiometry = constants(14);

free_BPTA = zeros(size(data(:,1)));
free_additive = zeros(size(data(:,1)));
Ke_BPTA = zeros(size(data(:,1)));
Kn_BPTA = zeros(size(data(:,1)));
Ke_additive = zeros(size(data(:,1)));
Kn_additive = zeros(size(data(:,1)));
sigma_BPTA = zeros(size(data(:,1)));
sigma_additive = zeros(size(data(:,1)));
Kseq = zeros(size(data(:,1)));

for i = 1:length(data(:,1))
    f_additive = data(i,1);
    T = data(i,2)+273.15;
    c_BPTA = (1-f_additive)*c_tot;
    c_additive = f_additive*c_tot;
    
    Ge_BPTA = H_BPTA - T*S_BPTA;
    Gn_BPTA = Ge_BPTA + NP_BPTA;
    Ge_additive = H_additive - T*S_additive;
    Gn_additive = Ge_additive + NP_additive;
    
    G_seq = H_seq - T*S_seq;
    
    Ke_BPTA(i) = exp(-Ge_BPTA./(R*T));
    Kn_BPTA(i) = exp(-Gn_BPTA./(R*T));
    Ke_additive(i) = exp(-Ge_additive./(R*T));
    Kn_additive(i) = exp(-Gn_additive./(R*T));
    sigma_BPTA(i) = Kn_BPTA(i)./Ke_BPTA(i);
    sigma_additive(i) = Kn_additive(i)./Ke_additive(i);
    Kseq(i) = exp(-G_seq./(R*T));

    [free_BPTA(i), free_additive(i), count] = calculate_mass_balance(Ke_BPTA(i), sigma_BPTA(i), Ke_additive(i), sigma_additive(i), Kseq(i), c_BPTA, c_additive, sequestrating_stoichiometry);
end
polymerized_BPTA = -sigma_BPTA.*free_BPTA + sigma_BPTA.*free_BPTA./(1-Ke_BPTA.*free_BPTA).^2;
polymerized_additive = -sigma_additive.*free_additive + sigma_additive.*free_additive./(1-Ke_additive.*free_additive).^2;

if sequestrating_stoichiometry == 1
    sequestrated_BPTA = free_BPTA.*free_additive.*Kseq;
    sequestrated_additive = free_BPTA.*free_additive.*Kseq;
elseif sequestrating_stoichiometry == 2
    sequestrated_BPTA = free_BPTA.*free_additive.*Kseq + 2.*free_BPTA.^2.*free_additive.*Kseq.^2;
    sequestrated_additive = free_BPTA.*free_additive.*Kseq + free_BPTA.^2.*free_additive.*Kseq.^2;
elseif sequestrating_stoichiometry == 3
    sequestrated_BPTA = free_BPTA.*free_additive.*Kseq + 2.*free_BPTA.^2.*free_additive.*Kseq.^2 + 3.*free_BPTA.^3.*free_additive.*Kseq.^3;
    sequestrated_additive = free_BPTA.*free_additive.*Kseq + free_BPTA.^2.*free_additive.*Kseq.^2 + free_BPTA.^3.*free_additive.*Kseq.^3;

end

CDcalc = (free_BPTA+free_additive).*Emono + polymerized_BPTA.*Epol_a + polymerized_additive.*Epol_b + Eseq.*(sequestrated_BPTA+sequestrated_additive);
end

