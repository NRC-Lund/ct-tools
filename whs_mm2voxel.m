function [ML, AP, DV] = whs_mm2voxel(ML_mm, AP_mm, DV_mm)
%WHS_MM2VOXEL   Convert WHS SD rat v1.01 voxels to millimeters.
%
% This function takes WHS SD rat v1.01 atlas voxel coordinates and converts
% them to millimeters (Papp et al, 2014 and 2015). The WHS SD rat v1.01 
% coordinate system is defined as Right, Anterior, Superior, with the mm
% origin at the decussation of the anterior commissure. Note that the voxel
% value is 1 less than the slice number.
%
% SYNTAX:
%   [ML, AP, DV] = whs_mm2voxel(ML_mm, AP_mm, DV_mm)
%
% INPUT:
%   ML_mm   - Column vector of coordinates in the medial-lateral direction.
%             Right is positive. Unit in millimeters.
%   AP_mm   - Column vector of coordinates in the anterior-posterior 
%             direction. Anterior is positive.  Unit in millimeters.
%   DV_mm   - Column vector of coordinates in the dorsal-ventral direction.
%             Dorsal is positive. Unit in millimeters.
%
% OUTPUT:
%   ML      - Column vector of coordinates in the medial-lateral direction.
%             Right is positive. Unit in voxels.
%   AP      - Column vector of coordinates in the anterior-posterior 
%             direction. Anterior is positive.  Unit in voxels.
%   DV      - Column vector of coordinates in the dorsal-ventral direction.
%             Dorsal is positive. Unit in voxels.
%

% Affine transformation, augmented matrix form: (Obtained from
% https://www.nitrc.org/docman/view.php/1081/2095/Coordinates_v1-v1.01.pdf)
tform = [0.0390625 0 0 -24.3359; ...
         0 0.0390625 0 -9.5312; ...
         0 0 0.0390625 -9.6875; ...
         0 0 0 1];
% Calculate:
co = tform \ [AP_mm(:) ML_mm(:) DV_mm(:) ones(numel(AP_mm),1)]';
ML = co(2,:)';
AP = co(1,:)';
DV = co(3,:)';
