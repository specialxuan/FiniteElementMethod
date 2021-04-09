function fePlot(x1, x2, x3, y1, y2, y3, var)
    mycolour = colour(var);
    for i = [1:length(var)]
        fill([x1(i), x2(i), x3(i)], [y1(i),y2(i), y3(i)], mycolour(i, :));
        hold on;
    end
end