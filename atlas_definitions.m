function defs = atlas_definitions()

% Paxinos & Watson: The rat brain 6th edition
ix = 1;
defs(ix).name = 'Paxinos and Watson: The Rat Brain (6th edition)';
defs(ix).dir = 'Paxinos_and_Watson_2007_The_Rat_Brain_6th_edition';
defs(ix).landmarks = 'atlas-landmarks.txt';
defs(ix).format = 'slices';
defs(ix).lambda_bregma = 8.7772; % rat
files = struct();
files(1).name = 'Paxinos_coronal_6th_edition.mat';
files(1).type = 'coronal';
files(2).name = 'Paxinos_sagittal_6th_edition.mat';
files(2).type = 'sagittal';
defs(ix).files = files;

% Waxholm Space atlas: SD rat v1.01
ix = 2;
defs(ix).name = 'Waxholm SD v1.01';
defs(ix).dir = 'WHS_SD_rat_v1.01';
defs(ix).landmarks = 'atlas-landmarks.txt';
defs(ix).format = 'volume';
defs(ix).lambda_bregma = 8.7772; % copied from Paxinos assuming rat strain makes no difference
files = struct();
files(1).name = 'WHS_SD_rat_T2star_v1.01.nii';
files(1).type = 'NIfTI';
files(2).name = 'segmented-atlas.mat';
files(2).type = 'segmented volume';
defs(ix).files = files;

% Franklin & Paxinos: The mouse brain 3rd edition
ix = 3;
defs(ix).name = 'Franklin and Paxinos: The Mouse Brain (3rd edition)';
defs(ix).dir = 'The mouse brain in stereotaxic coordinates Third edition Paxinos';
defs(ix).landmarks = 'atlas-landmarks.txt';
defs(ix).format = 'slices';
defs(ix).lambda_bregma = 4.21; % mouse
files = struct();
files(1).name = 'Mouse_Paxinos_Coronal_3rd_Edition_150dpi.mat';
files(1).type = 'coronal';
defs(ix).files = files;
