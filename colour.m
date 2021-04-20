function mycolour = colour(var)
    max_var = max(abs(var));
    green = var / max_var;
    green(green < 0) = 0;
    red = var / max_var;
    red(red > 0) = 0;
    red = - red;
    mycolour = [red, green, 1-green - red];
end
