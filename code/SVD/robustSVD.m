%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage:
% matlab robustSVD.m
% Function:
% perform the gene enrichment analysis by the Wilcoxon scoring method
% on the polarized gene eigenvectors.
% Input:
% - "data/CREF_Matrix/{species}_CREF_Matrix.tsv
% - "data/geneName/{species}_genename.txt"
% - "data/motif_symbol.txt"
% Output:
% - "results/SVD/{species}_S.txt"
% - "results/SVD/{species}_U0_9.txt"
% - "results/SVD/{species}_V0_1402.txt"
% - "results/loadings/{species}_{gene_*}_loadings.txt"
% - "results/loadings/{species}_{motif_*}_loadings.txt"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd SVD
% pwd = './code/SVD'

% add MATLAB search path
addpath(genpath('.'))

config = yaml.ReadYaml('../config.yaml');

%% directory settings
data_dir = '../../data';
result_dir = '../../results';
motif_path = fullfile(data_dir, 'motif_symbol.txt');
CisMatrix_dir = fullfile(data_dir, 'CREF_Matrix');
geneName_dir = fullfile(data_dir, 'geneName');
SVD_dir = fullfile(result_dir, 'SVD');
if ~exist(SVD_dir, 'dir')
    mkdir(SVD_dir);
end
loading_dir = fullfile(result_dir, 'loadings', filesep);
if ~exist(loading_dir, 'dir')
    mkdir(loading_dir);
end

%% parameters settings
speciesNum = length(config.species_list);
species_list = {};
for i = 1:speciesNum
    species_list{i} = config.species_list{i};
end

%% run robust SVD
UCell = cell(1, speciesNum);
VCell = cell(1, speciesNum);
SCell = cell(1, speciesNum);
LCell = cell(1, speciesNum);
sparseCell = cell(1, speciesNum);
LCell_fusion = cell(1, speciesNum);

for i = 1:speciesNum
    tic;
    disp(species_list{i});
    CisMatrix = load(fullfile(CisMatrix_dir, [species_list{i}, '_CREF_Matrix.tsv']));
    [UCell{i}, SCell{i}, VCell{i}, LCell{i}, sparseCell{i}] = robustPCA(CisMatrix);
    toc;
end

%% correct sign
for j = 1:1403
    if min(VCell{1}(:, j)) + max(VCell{1}(:, j)) < 0
        VCell{1}(:, j) = -VCell{1}(:, j);
        UCell{1}(:, j) = -UCell{1}(:, j);
    end
    for i = 2:speciesNum
        if(dot(VCell{1}(:, j), VCell{i}(:, j)) < 0)
            VCell{i}(:, j) = -VCell{i}(:, j);
            UCell{i}(:, j) = -UCell{i}(:, j);
        end
    end
end

%% save UCell, VCell ,SCell of each species to .txt files separately
for i = 1:speciesNum
    filename = fullfile(SVD_dir, [species_list{i}, '_U0_9.txt']);
    cell2txt(UCell{i}, filename, 10);

    filename = fullfile(SVD_dir, [species_list{i}, '_V0_1402.txt']);
    cell2txt(VCell{i}, filename, 1403);

    filename = fullfile(SVD_dir, [species_list{i}, '_S.txt']);
    cell2txt(diag(SCell{i}), filename);
end

% save the matrix of each species into a .mat file alone
% (requires a large amount of storage space).
% for i = 1:speciesNum
%     U = UCell{i};
%     % save(fullfile(SVD_dir, [species_list{i}, '_U.mat']), 'U', '-v7.3');
%     save(fullfile(SVD_dir, [species_list{i}, '_U.mat']), 'U');
%     S = SCell{i};
%     save(fullfile(SVD_dir, [species_list{i}, '_S.mat']), 'S');
%     V = VCell{i};
%     save(fullfile(SVD_dir, [species_list{i}, '_V.mat']), 'V');
%     L = LCell{i};
%     save(fullfile(SVD_dir, [species_list{i}, '_L.mat']), 'L');
%     sparse = sparseCell{i};
%     save(fullfile(SVD_dir, [species_list{i}, '_sparse.mat']), 'sparse');
% end

%% save motif and gene loadings
geneWrite(species_list, UCell, 10, loading_dir, geneName_dir);
motifWrite(species_list, VCell, 10, loading_dir, motif_path);
