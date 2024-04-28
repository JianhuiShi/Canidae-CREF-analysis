function cell2txt(onecell, filename, varargin)

fid = fopen(filename,'w');
[nrow, ncol] = size(onecell);

if length(varargin) == 1
    ncol = varargin{1};
end

for j = 1:nrow
    for k = 1:ncol-1
        fprintf(fid,'%.12f\t',onecell(j,k));
    end
    fprintf(fid, '%.12f\n', onecell(j, ncol));
end
fclose(fid);
