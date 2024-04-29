% Copyright © 2024, Jianhui Shi & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

function geneWrite(species, UCell, col, loading_dir, geneName_dir)

for i = 1:length(species)
    genename = readtable(fullfile(geneName_dir, [species{i}, '_genename.txt']), 'ReadVariableNames', false);
    genename = genename.Var1;
    for k = 1:col
        filename = fullfile(loading_dir, [species{i}, '_gene_', num2str(k - 1), '_loadings.txt']);
        fid = fopen(filename,'w');
        column = UCell{i}(:, k);
        [~, y] = sort(column, 'descend');
        fprintf(fid, '%s\t', 'Symbol');
        fprintf(fid, '%s\n', 'Loading');
        % 'column' are printed in descending order to the file
        for j = 1:length(column)
            fprintf(fid, '%s\t', genename{y(j)});
            fprintf(fid, '%f\n', column(y(j)));
        end
        fclose(fid);
    end
end
end