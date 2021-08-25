function cost = calculate_cost_function(parameters, data, constants)

[CDcalc] = simulate_polymerization_sequestration(parameters,constants,data);
cost = (data(:,3)-CDcalc);

end