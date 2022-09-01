function [ML_mm, AP_mm, DV_mm] = whs_voxel2mm(ML, AP, DV)
%WHS_VOXEL2MM   Convert WHS SD rat v1.01 voxels to millimeters.
%
% This function takes WHS SD rat v1.01 atlas voxel coordinates and converts
% them to millimeters (Papp et al, 2014 and 2015). The WHS SD rat v1.01 
% coordinate system is defined as Right, Anterior, Superior, with the mm
% origin at the decussation of the anterior commissure. Note that the voxel
% value is 1 less than the slice number.
%
% SYNTAX:
%   [ML_mm, AP_mm, DV_mm] = whs_voxel2mm(ML, AP, DV)
%   tform = whs_voxel2mm()
%
% INPUT:
%   ML      - Column vector of coordinates in the medial-lateral direction.
%             Right is positive. Unit in voxels.
%   AP      - Column vector of coordinates in the anterior-posterior 
%             direction. Anterior is positive.  Unit in voxels.
%   DV      - Column vector of coordinates in the dorsal-ventral direction.
%             Dorsal is positive. Unit in voxels.
%
% OUTPUT:
%   ML_mm   - Column vector of coordinates in the medial-lateral direction.
%             Right is positive. Unit in millimeters.
%   AP_mm   - Column vector of coordinates in the anterior-posterior 
%             direction. Anterior is positive.  Unit in millimeters.
%   DV_mm   - Column vector of coordinates in the dorsal-ventral direction.
%             Dorsal is positive. Unit in millimeters.
%   tform   - Affine transformation.

% Affine transformation, augmented matrix form: (Obtained from
% https://www.nitrc.org/docman/view.php/1081/2095/Coordinates_v1-v1.01.pdf)
tform = [0.0390625 0 0 -24.3359; ...
         0 0.0390625 0 -9.5312; ...
         0 0 0.0390625 -9.6875; ...
         0 0 0 1];
if nargin==0
    ML_mm = tform;
    ML_mm(1,4) = tform(2,4);
    ML_mm(2,4) = tform(1,4);
    return
end
% Calculate:
co = tform * [AP(:) ML(:) DV(:) ones(numel(AP),1)]';
ML_mm = co(2,:)';
AP_mm = co(1,:)';
DV_mm = co(3,:)';
