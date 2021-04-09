function mycolour = colour(var)
    max_var = max(abs(var));
    red = abs(var / max_var);
    mycolour = [red, zeros(size(red)), 1-red];
end
