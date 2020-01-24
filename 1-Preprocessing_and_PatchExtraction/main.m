% Step 0 : [Assumed to be completed]
% TCGA data should be downloaded and stored in the following file struct:
%   /path/to/imgFolder
%       ../TCGA ID
%           ../ contains all svs and xml files
%
% for example:
%   /path/to/final_cohort
%       ../TCGA-2f-A9KO
%           contains:
%           TCGA-2F-A9KO-01Z-00-DX1.195576CF-B739-4BD9-B15B-4A70AE287D3E.svs
%           TCGA-2F-A9KO-01Z-00-DX1.195576CF-B739-4BD9-B15B-4A70AE287D3E.xml

% Step 1: pre-process images and extract patches 
addpath(genpath('/path/to/dependenies'))
imgDir  = '/path/to/imgFolder';
patchDir = '/path/where/to/save/patches';
svsList = dir([imgDir '/**/*.svs']);
outcomeFile = '/path/to/outcomes/file';
refimg = load('/path/to/ref/TCGA-4Z-AA81-01Z-00-DX1_ref.mat');

for i = 1:numel(svsList)
    xml_file = [svsList(i).folder filesep strrep(svsList(i).name,'svs','xml')];
    svs_file = [svsList(i).folder filesep svsList(i).name];
    tcga_id = strsplit(strrep(svsList(i).name,'.svs',''),'_');
    tcga_id = strjoin(tcga_id(1:3),'_');
    savename = strrep(svsList(i).name,'.svs','');
    load_disp_tcga(xml_file, svs_file, imgDir, tcga_id, savename)
    tcga_newnorm(refimg,[imgDir filesep tcga_id filesep svs_file])
    TCGA_pullBoxes(strrep([imgDir filesep tcga_id filesep svs_file],'.svs',''),imgDir,patchDir,outcomeFile)
    imglist = dir(['/data/MIP/TCGABLCA/classification/gtom_FINAL/prelim' filesep tcga_id filesep 'decon' filesep '*.jpeg']);
    
end