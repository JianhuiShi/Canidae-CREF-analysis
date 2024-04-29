% Copyright © 2024, Jianhui Shi & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

function motifWrite(species, VCell, col, loading_dir, motif_path)
motif = readtable(motif_path, 'ReadVariableNames', false, 'Delimiter', ' ');
motif = motif.Var1;
motif{1404} = '';

for i = 1:length(species)
    for k = 1:col
        filename = fullfile(loading_dir, [species{i}, '_motif_', num2str(k - 1), '_loadings.txt']);
        fid = fopen(filename,'w');
        column = VCell{i}(:,k);
        [~, y] = sort(column, 'descend');
        fprintf(fid, '%s\t', 'Symbol');
        fprintf(fid, '%s\n', 'Loading');
        % 'column' are printed in descending order to the file
        for j = 1:length(column)
            fprintf(fid, '%s\t', motif{y(j)});
            fprintf(fid, '%f\n', column(y(j)));
        end
        fclose(fid);
    end
end
end