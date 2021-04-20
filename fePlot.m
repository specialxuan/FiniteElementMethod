function fePlot(result, var)
    mycolour = colour(var);
    for i = [1:length(var)]
        fill([result.x1(i), result.x2(i), result.x3(i)], [result.y1(i), result.y2(i), result.y3(i)], mycolour(i, :), 'Edgecolor', 'none');
        hold on;
    end
    axis square;
end
