function structure(m, n)
    inputFile = fopen("fe.csv")
    fprintf(inputFile, 'INPUT, Degree of freedom is %d', 2 * m * n)
end
