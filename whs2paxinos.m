function [ML_pax, AP_pax, DV_pax] = whs2paxinos(ML_whs, AP_whs, DV_whs, type)
%WHS2PAXINOS   Convert WHS SD rat v1.01 coordinates to Paxinos coordinates.
%
% This function takes coordinates given in mm in the WHS SD rat v1.01 atlas
% coordinate system (Papp et al, 2014 and 2015) and converts them to the 
% Paxinos coordinate system (Paxinos & Watson, 2007). The WHS SD rat v1.01 
% coordinate system is defined as Right, Anterior, Superior, with the 
% origin at the decussation of the anterior commissure. Units must be in 
% millimeters. The Paxinos system is also defined as RAS, but with the
% origin at bregma.
%
% There is not one optimal way to do this transformation, and several
% transformation types can be selected with the 'type' argument. See
% get_whs2paxinos_tform for more details.
%
% SYNTAX:
%   [ML_pax, AP_pax, DV_pax] = whs2paxinos(ML_whs, AP_whs, DV_whs)
%   [ML_pax, AP_pax, DV_pax] = whs2paxinos(ML_whs, AP_whs, DV_whs, type)
%
% INPUT:
%   ML_whs  - Column vector of coordinates in the medial-lateral direction.
%             Right is positive. Unit in millimeters.
%   AP_whs  - Column vector of coordinates in the anterior-posterior 
%             direction. Anterior is positive.  Unit in millimeters.
%   DV_whs  - Column vector of coordinates in the dorsal-ventral direction.
%             Dorsal is positive. Unit in millimeters.
%   type    - (optional, default='Anisotropic1') See get_whs2paxinos_tform 
%             for more details.
%
% OUTPUT:
%   ML_pax  - Column vector of coordinates in the medial-lateral direction.
%             Right is positive. Unit in millimeters.
%   AP_pax  - Column vector of coordinates in the anterior-posterior 
%             direction. Anterior is positive.  Unit in millimeters.
%   DV_pax  - Column vector of coordinates in the dorsal-ventral direction.
%             Dorsal is positive. Unit in millimeters.

% Get affine transformation:
if nargin==3
    type = 'Anisotropic1';
end
tform = get_whs2paxinos_tform(type); % Augmented matrix form.

% Calculate:
co = tform * [ML_whs(:) AP_whs(:) DV_whs(:) ones(numel(ML_whs),1)]';
ML_pax = co(1,:)';
AP_pax = co(2,:)';
DV_pax = co(3,:)';

% Reshape:
ML_pax = reshape(ML_pax,size(ML_whs));
AP_pax = reshape(AP_pax,size(AP_whs));
DV_pax = reshape(DV_pax,size(DV_whs));