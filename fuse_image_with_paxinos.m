function [ImgCoronal, RA_Coronal, coronal, ImgSagittal, RA_Sagittal, sagittal] = ...
    fuse_image_with_paxinos(Img, tform, DV_step, ML_step, Animal)
%FUSE_IMAGE_WITH_PAXINOS   Fuse image with Paxinos atlas.
%
% Creates a fused Image from an input image and an affine transformation.
%
% SYNTAX:
%   [ImgCoronal, RA_Coronal, coronal] = ...
%           fuse_image_with_paxinos(img, tform)
%   [ImgCoronal, RA_Coronal, coronal, ImgSagittal, RA_Sagittal, sagittal] = ...
%           fuse_image_with_paxinos(img, tform)
%   [...] = fuse_image_with_paxinos(img, tform, DV_step)
%   [...] = fuse_image_with_paxinos(img, tform, DV_step, ML_step)
%   [...] = fuse_image_with_paxinos(img, tform, DV_step, ML_step, Animal)
%
% INPUT:
%   img         - Input image [ML,AP,DV] with increasing indices in the
%                 directions right, anterior and dorsal/superior.
%   tform       - Affine transformation on augmented matrix form or an
%                 affine3d object.
%   DV_step     - (optional, default=1) DV step size for coronal view.
%   ML_step     - (optional, default=1) ML step size for sagittal view.
%   Animal      - (optional, default='rat') animal species (to select atlas), 'rat' or 'mouse' 
%
% OUTPUT:
%   ImgCoronal  - Fused coronal image [ML,DV,RGB,AP].
%   RA_Coronal  - Reference object for coronal images.
%   coronal     - Structure with coronal Paxinos images.
%   ImgSagittal - Fused sagittal image [AP,DV,RGB,ML].
%   RA_Sagittal - Reference object for sagittal images.
%   sagittal    - Structure with sagittal Paxinos images.
%

% Step size:
if nargin<3
    DV_step = 1;
end
if nargin<4
    ML_step = 1;
end
if nargin<5
    Animal = 'rat';
end

% Load Paxinos:
if Animal == "rat"
  fprintf('Loading Paxinos Rat atlas...')
  dname = fullfile(fileparts(mfilename('fullpath')), 'Atlases', ...
      'Paxinos_and_Watson_2007_The_Rat_Brain_6th_edition'); % Get Rat Paxinos directory.
  coronal = load(fullfile(dname,'Paxinos_coronal_6th_edition.mat'));
  sagittal = load(fullfile(dname,'Paxinos_sagittal_6th_edition.mat'));
  % Limits of output images in Paxinos mm:
  MlLim = [-9 9];
  ApLim = [-20 10];
  DvLim = [-15 2]; % for rat!
elseif Animal == "mouse"
  fprintf('Loading Paxinos Mouse atlas...')
  dname = fullfile(fileparts(mfilename('fullpath')), 'Atlases', ...
      'The mouse brain in stereotaxic coordinates Third edition Paxinos'); % Get Mouse Paxinos directory.
  coronal = load(fullfile(dname,'Mouse_Paxinos_Coronal_3rd_Edition_150dpi.mat'));
  % Limits of output images in Paxinos mm:
  MlLim = [-6 6];
  ApLim = [-8.5 4.5];
  DvLim = [-6.2 1]; % for mouse
else
  error('Unknown animal type!');
end
fprintf('done.\n')

% Create affine3d object:
if ~isa(tform, 'affine3d')
    tform = affine3d(tform');
end

% Decompose transform (we need to extract the scaling):
[tform_scale, tform_rot, tform_trans] = decompose_affine_tform(tform.T');
tform_rigid = tform_rot*tform_trans; % If we do scaling to mm with imwarp we will reduce the resolution, so we do scaling by changing the world reference frame instead.
tform_rigid = affine3d(tform_rigid');

% Warp:
fprintf('Warping image...')
ImgReg = permute(Img, [2 1 3]); % Need to have [AP,ML,DV], because intrinsic coordinates relates to image indices as [y,x,z].
RA_3d = imref3d(size(ImgReg));
[ImgReg,RA_3d_Reg_Unscaled] = imwarp(ImgReg, RA_3d, tform_rigid); % In the output ref object, the world coordinates are the old subscripts.
fprintf('done.\n')

% Creating reference object with scaled world coordinates:
RA_3d_Reg = RA_3d_Reg_Unscaled;
[x,y,z] = tform.transformPointsInverse(0,0,0); % Find new origo in old coordinates.
[x,y,z] = tform_rigid.transformPointsForward(x,y,z); % Transform to new unscaled space.
RA_3d_Reg.XWorldLimits = RA_3d_Reg.XWorldLimits - x; % Translate so that origo is zero before we scale...
RA_3d_Reg.YWorldLimits = RA_3d_Reg.YWorldLimits - y;
RA_3d_Reg.ZWorldLimits = RA_3d_Reg.ZWorldLimits - z;
RA_3d_Reg.XWorldLimits = tform_scale(1)*RA_3d_Reg.XWorldLimits; % Scale...
RA_3d_Reg.YWorldLimits = tform_scale(1)*RA_3d_Reg.YWorldLimits;
RA_3d_Reg.ZWorldLimits = tform_scale(1)*RA_3d_Reg.ZWorldLimits;

% Fuse coronal:
disp('Fusing coronal images with atlas...');
SliceIx = 1:DV_step:height(coronal.T);
iSlice = round(numel(SliceIx)/2); % Take a middle slice to estimate the size of the final image.
[AP_Index, ML_Index, DV_Index] = ...
    RA_3d_Reg.worldToSubscript(0,coronal.T.AP(SliceIx(iSlice)),0);  % Get AP slice number.
ImgReg2 = squeeze(ImgReg(AP_Index,:,:));                            % Native AP slice.
ImgReg2 = permute(ImgReg2, [2 1]);                                  % Now in [DV,ML]
RA_2d = imref2d(size(ImgReg2), RA_3d_Reg.XWorldLimits, RA_3d_Reg.ZWorldLimits); % Reference object.
ImgPax = coronal.Img(:,:,SliceIx(iSlice));      % Paxinos image.
[ImgPax, RA_fused] = imfuse(ImgPax, coronal.RA(SliceIx(iSlice)), ...
    ImgReg2, RA_2d, 'Scaling', 'none');   % Fuse.
[x,y] = RA_fused.worldToIntrinsic(MlLim, DvLim);    % Intrinsic coordinates.
ImgSize = [round(diff(y)) round(diff(x))];          % Now we have the output size.
ImgCoronal = repmat(ImgPax(1),ImgSize(1),ImgSize(2),3,numel(SliceIx));   % Allocate variable.
for iSlice = 1:numel(SliceIx)                       % Do all slices...
    fprintf('   Slice %d/%d\n', iSlice, numel(SliceIx))
    [AP_Index, ML_Index, DV_Index] = ...
        RA_3d_Reg.worldToSubscript(0,coronal.T.AP(SliceIx(iSlice)),0);  % Get AP slice number.
    ImgReg2 = squeeze(ImgReg(AP_Index,:,:));    % Native AP slice.
    ImgReg2 = permute(ImgReg2, [2 1]);          % Now in [DV,ML]
    RA_2d = imref2d(size(ImgReg2), RA_3d_Reg.XWorldLimits, RA_3d_Reg.ZWorldLimits); % Reference object.
    ImgPax = coronal.Img(:,:,SliceIx(iSlice));  % Paxinos image.
    ImgPax = imcomplement(ImgPax);              % Dark background.
    ImgPax = imadjust(ImgPax);                  % Adjust contrast.
    [ImgPax, RA_fused] = imfuse(ImgPax, coronal.RA(SliceIx(iSlice)), ...
        ImgReg2, RA_2d, 'Scaling', 'none');     % Fuse.
    OV = affineOutputView([size(ImgPax,1) size(ImgPax,2)], affine2d); % Define OutputView.
    OV.ImageSize = ImgSize;
    OV.XWorldLimits = MlLim;
    OV.YWorldLimits = DvLim;
    [ImgCoronal(:,:,:,iSlice),RA_Coronal] = ...
        imwarp(ImgPax, RA_fused, affine2d, 'OutputView', OV);
end

% Fuse sagittal:
if nargout>3
    disp('Fusing sagittal images with atlas...');
    SliceIx = 1:ML_step:height(sagittal.T);
    iSlice = 1;
    [AP_Index, ML_Index, DV_Index] = ...
        RA_3d_Reg.worldToSubscript(sagittal.T.ML(SliceIx(iSlice)),0,0); % Get ML slice number.
    ImgReg2 = squeeze(ImgReg(:,ML_Index,:)); % Native ML slice. Was [AP,ML,DV], now [AP,DV].
    ImgReg2 = permute(ImgReg2, [2 1]);       % Now in [DV,AP]
    RA_2d = imref2d(size(ImgReg2), RA_3d_Reg.YWorldLimits, RA_3d_Reg.ZWorldLimits); % Reference object.
    ImgPax = sagittal.Img(:,:,SliceIx(iSlice));     % Paxinos image.
    [ImgPax, RA_fused] = imfuse(ImgPax, sagittal.RA(SliceIx(iSlice)), ...
        ImgReg2, RA_2d, 'Scaling', 'none'); % Fuse.
    [x,y] = RA_fused.worldToIntrinsic(ApLim, DvLim);    % Intrinsic coordinates.
    ImgSize = [round(diff(y)) round(diff(x))];          % Now we have the output size.
    ImgSagittal = repmat(ImgPax(1),ImgSize(1),ImgSize(2),3,numel(SliceIx));   % Allocate variable.
    for iSlice = 1:numel(SliceIx)                       % Do all slices...
        fprintf('   Slice %d/%d\n', iSlice, numel(SliceIx))
        [AP_Index, ML_Index, DV_Index] = ...
            RA_3d_Reg.worldToSubscript(sagittal.T.ML(SliceIx(iSlice)),0,0);
        ImgReg2 = squeeze(ImgReg(:,ML_Index,:));
        ImgReg2 = permute(ImgReg2, [2 1]);  % Now in [DV,ML]
        RA_2d = imref2d(size(ImgReg2), RA_3d_Reg.YWorldLimits, RA_3d_Reg.ZWorldLimits);
        ImgPax = sagittal.Img(:,:,SliceIx(iSlice));
        ImgPax = imcomplement(ImgPax);
        ImgPax = imadjust(ImgPax);
        [ImgPax, RA_fused] = imfuse(ImgPax, sagittal.RA(SliceIx(iSlice)), ...
            imadjust(ImgReg2), RA_2d, 'Scaling', 'none');
        OV = affineOutputView([size(ImgPax,1) size(ImgPax,2)], affine2d); % Define OutputView.
        OV.ImageSize = ImgSize;
        OV.XWorldLimits = ApLim;
        OV.YWorldLimits = DvLim;
        [ImgSagittal(:,:,:,iSlice),RA_Sagittal] = ...
            imwarp(ImgPax, RA_fused, affine2d, 'OutputView', OV);
    end
end
