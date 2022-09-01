function [coronal_out, sagittal_out] = fuse_image_with_atlas(Img, tform, varargin)
%FUSE_IMAGE_WITH_PAXINOS   Fuse image with Paxinos atlas.
%
% Creates a fused Image from an input image and an affine transformation.
%
% SYNTAX:
%   [ImgCoronal, RA_Coronal, coronal] = ...
%           fuse_image_with_paxinos(img, tform)
%   [ImgCoronal, RA_Coronal, coronal, ImgSagittal, RA_Sagittal, sagittal] = ...
%           fuse_image_with_paxinos(img, tform)
%   [...] = fuse_image_with_paxinos(img, tform, AP_step)
%   [...] = fuse_image_with_paxinos(img, tform, AP_step, ML_step)
%
% INPUT:
%   img         - Input image [ML,AP,DV] with increasing indices in the
%                 directions right, anterior and dorsal/superior.
%   tform       - Affine transformation on augmented matrix form or an
%                 affine3d object.
%   AP_step     - (optional, default=1) AP step size for coronal view.
%   ML_step     - (optional, default=1) ML step size for sagittal view.
%
% OUTPUT:
%   ImgCoronal  - Fused coronal image [ML,DV,RGB,AP].
%   RA_Coronal  - Reference object for coronal images.
%   coronal     - Structure with coronal Paxinos images.
%   ImgSagittal - Fused sagittal image [AP,DV,RGB,ML].
%   RA_Sagittal - Reference object for sagittal images.
%   sagittal    - Structure with sagittal Paxinos images.
%

% Optional parameters:
coronal = [];
sagittal = [];
AP_step = 1;        % Number of coronal slices to skip.
ML_step = 1;        % Number of sagittal slices to skip.
ML_lim = [-9 9];    % Limits of output images in atlas coordinates...
AP_lim = [-20 10];  % ...
DV_lim = [-15 2];   % ...
if mod(length(varargin),2)
    error('Incomplete property-value pairs!');
else
    % For every pair of input arguments ...
    for i = 1:2:length(varargin)
        % Check which property is to be set
        switch lower(varargin{i})
            case 'coronal'
                coronal = varargin{i+1};
            case 'sagittal'
                sagittal = varargin{i+1};
            case 'apstep'
                AP_step = varargin{i+1};
            case 'mlstep'
                ML_step = varargin{i+1};
            case 'mllim'
                ML_lim = varargin{i+1};
            case 'aplim'
                AP_lim = varargin{i+1};
            case 'dvlim'
                DV_lim = varargin{i+1};
            otherwise % Invalid property
                error(['The property ',varargin{i},' is invalid!']);
        end
    end
end
if isempty(coronal) && isempty(sagittal)
    error('Both coronal and sagittal slices are empty.')
end

% % Load Paxinos:
% fprintf('Loading Paxinos atlas...')
% dname = fullfile(fileparts(mfilename('fullpath')), 'Atlases', ...
%     'Paxinos_and_Watson_2007_The_Rat_Brain_6th_edition'); % Get Paxinos directory.
% coronal = load(fullfile(dname,'Paxinos_coronal_6th_edition.mat'));
% sagittal = load(fullfile(dname,'Paxinos_sagittal_6th_edition.mat'));
% fprintf('done.\n')

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
if isempty(coronal)
    coronal_out = [];
else
    disp('Fusing coronal images with atlas...');
    SliceIx = 1:AP_step:height(coronal.T);
    iSlice = 1;
    fprintf('   Slice %d/%d\n', iSlice, numel(SliceIx))
    [AP_Index, ML_Index, DV_Index] = ...
        RA_3d_Reg.worldToSubscript(0,coronal.T.AP(SliceIx(iSlice)),0);  % Get AP slice number.
    ImgReg2 = squeeze(ImgReg(AP_Index,:,:));                            % Native AP slice.
    ImgReg2 = permute(ImgReg2, [2 1]);                                  % Now in [DV,ML]
    RA_2d = imref2d(size(ImgReg2), RA_3d_Reg.XWorldLimits, RA_3d_Reg.ZWorldLimits); % Reference object.
    ImgPax = coronal.Img(:,:,SliceIx(iSlice));      % Paxinos image.
    ImgPax = imcomplement(ImgPax);                  % Dark background.
    ImgPax = imadjust(ImgPax);                      % Adjust contrast.
    [ImgPax, RA_fused] = imfuse(ImgPax, coronal.RA(SliceIx(iSlice)), ...
        imadjust(ImgReg2), RA_2d, 'Scaling', 'none'); % Fuse.
    [ix1,ix2] = RA_fused.worldToSubscript(ML_lim, DV_lim);    % Crop...
    if isnan(ix1(1)); ix1(1)=1;end                          % ...
    if isnan(ix2(1)); ix2(1)=1;end                          % ...
    if isnan(ix1(2)); ix1(2)=RA_fused.ImageSize(1);end      % ...
    if isnan(ix2(2)); ix2(2)=RA_fused.ImageSize(2);end      % ...
    ImgPax = ImgPax(ix1(1):ix1(2),ix2(1):ix2(2),:);         % ...
    ImgSize = size(ImgPax);
    ImgOut = repmat(ImgPax,1,1,1,numel(SliceIx));
    for iSlice = 2:numel(SliceIx)
        fprintf('   Slice %d/%d\n', iSlice, numel(SliceIx))
        [AP_Index, ML_Index, DV_Index] = ...
            RA_3d_Reg.worldToSubscript(0,coronal.T.AP(SliceIx(iSlice)),0);
        ImgReg2 = squeeze(ImgReg(AP_Index,:,:));
        ImgReg2 = permute(ImgReg2, [2 1]);  % Now in [DV,ML]
        RA_2d = imref2d(size(ImgReg2), RA_3d_Reg.XWorldLimits, RA_3d_Reg.ZWorldLimits);
        ImgPax = coronal.Img(:,:,SliceIx(iSlice));
        ImgPax = imcomplement(ImgPax);
        ImgPax = imadjust(ImgPax);
        [ImgPax, RA_fused] = imfuse(ImgPax, coronal.RA(SliceIx(iSlice)), ...
            ImgReg2, RA_2d, 'Scaling', 'none');
        [ix1,ix2] = RA_fused.worldToSubscript(ML_lim, DV_lim);    % Crop...
        if isnan(ix1(1)); ix1(1)=1;end                          % ...
        if isnan(ix2(1)); ix2(1)=1;end                          % ...
        if isnan(ix1(2)); ix1(2)=RA_fused.ImageSize(1);end      % ...
        if isnan(ix2(2)); ix2(2)=RA_fused.ImageSize(2);end      % ...
        ImgPax = ImgPax(ix1(1):ix1(2),ix2(1):ix2(2),:);
        ImgOut(:,:,:,iSlice) = imresize(ImgPax, [ImgSize(1) ImgSize(2)]);
    end
    RA_Coronal = imref2d(ImgSize, ML_lim, DV_lim); % Reference object.
    coronal_out = struct('Img', ImgOut, 'RA', RA_Coronal, 'T', coronal.T);
end

% Fuse sagittal:
if isempty(sagittal)
    sagittal_out = [];
else
    disp('Fusing sagittal images with atlas...');
    SliceIx = 1:ML_step:height(sagittal.T);
    iSlice = 1;
    fprintf('   Slice %d/%d\n', iSlice, numel(SliceIx))
    [AP_Index, ML_Index, DV_Index] = ...
        RA_3d_Reg.worldToSubscript(sagittal.T.ML(SliceIx(iSlice)),0,0); % Get ML slice number.
    ImgReg2 = squeeze(ImgReg(:,ML_Index,:)); % Native ML slice. Was [AP,ML,DV], now [AP,DV].
    ImgReg2 = permute(ImgReg2, [2 1]);       % Now in [DV,AP]
    RA_2d = imref2d(size(ImgReg2), RA_3d_Reg.YWorldLimits, RA_3d_Reg.ZWorldLimits); % Reference object.
    ImgPax = sagittal.Img(:,:,SliceIx(iSlice));     % Paxinos image.
    ImgPax = imcomplement(ImgPax);                  % Dark background.
    ImgPax = imadjust(ImgPax);                      % Adjust contrast.
    [ImgPax, RA_fused] = imfuse(ImgPax, sagittal.RA(SliceIx(iSlice)), imadjust(ImgReg2), RA_2d); % Fuse.
    [ix1,ix2] = RA_fused.worldToSubscript(AP_lim, DV_lim);    % Crop...
    if isnan(ix1(1)); ix1(1)=1;end                          % ...
    if isnan(ix2(1)); ix2(1)=1;end                          % ...
    if isnan(ix1(2)); ix1(2)=RA_fused.ImageSize(1);end      % ...
    if isnan(ix2(2)); ix2(2)=RA_fused.ImageSize(2);end      % ...
    ImgPax = ImgPax(ix1(1):ix1(2),ix2(1):ix2(2),:);         % ...
    ImgSize = size(ImgPax);
    ImgOut = repmat(ImgPax,1,1,1,numel(SliceIx));
    for iSlice = 2:numel(SliceIx)
        fprintf('   Slice %d/%d\n', iSlice, numel(SliceIx))
        [AP_Index, ML_Index, DV_Index] = ...
            RA_3d_Reg.worldToSubscript(sagittal.T.ML(SliceIx(iSlice)),0,0);
        ImgReg2 = squeeze(ImgReg(:,ML_Index,:));
        ImgReg2 = permute(ImgReg2, [2 1]);  % Now in [DV,ML]
        RA_2d = imref2d(size(ImgReg2), RA_3d_Reg.YWorldLimits, RA_3d_Reg.ZWorldLimits);
        ImgPax = sagittal.Img(:,:,SliceIx(iSlice));
        ImgPax = imcomplement(ImgPax);
        ImgPax = imadjust(ImgPax);
        [ImgPax, RA_fused] = imfuse(ImgPax, sagittal.RA(SliceIx(iSlice)), imadjust(ImgReg2), RA_2d);
        [ix1,ix2] = RA_fused.worldToSubscript(AP_lim, DV_lim);    % Crop...
        if isnan(ix1(1)); ix1(1)=1;end                          % ...
        if isnan(ix2(1)); ix2(1)=1;end                          % ...
        if isnan(ix1(2)); ix1(2)=RA_fused.ImageSize(1);end      % ...
        if isnan(ix2(2)); ix2(2)=RA_fused.ImageSize(2);end      % ...
        ImgPax = ImgPax(ix1(1):ix1(2),ix2(1):ix2(2),:);         % ...
        ImgOut(:,:,:,iSlice) = imresize(ImgPax, [ImgSize(1) ImgSize(2)]);
    end
    RA_Sagittal = imref2d(ImgSize, AP_lim, DV_lim); % Reference object.
    sagittal_out = struct('Img', ImgOut, 'RA', RA_Sagittal, 'T', sagittal.T);
end


