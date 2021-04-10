function mycolour = colour(var)
    max_var = max(abs(var));
    red = var / max_var;
    red(red < 0) = 0;
    green = var / max_var;
    green(green > 0) = 0;
    green = - green;
    mycolour = [red, green, 1-red - green];
end
