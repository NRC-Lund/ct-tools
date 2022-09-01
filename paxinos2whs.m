function [ML_whs, AP_whs, DV_whs] = paxinos2whs(ML_pax, AP_pax, DV_pax, type)
%PAXINOS2WHS   Convert Paxinos coordinates to WHS SD rat v1.01 coordinates.
%
% This function takes coordinates given in mm in the Paxinos coordinate 
% system (bregma centered) and converts them to the WHS SD rat v1.01 atlas
% coordinate system (Papp et al, 2014 and 2015). The WHS SD rat v1.01 
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
%   [ML_whs, AP_whs, DV_whs] = paxinos2whs(ML_pax, AP_pax, DV_pax)
%   [ML_whs, AP_whs, DV_whs] = paxinos2whs(ML_pax, AP_pax, DV_pax, type)
%
% INPUT:
%   ML_pax  - Column vector of coordinates in the medial-lateral direction.
%             Right is positive. Unit in millimeters.
%   AP_pax  - Column vector of coordinates in the anterior-posterior 
%             direction. Anterior is positive.  Unit in millimeters.
%   DV_pax  - Column vector of coordinates in the dorsal-ventral direction.
%             Dorsal is positive. Unit in millimeters.
%   type    - (optional, default='Anisotropic1') See get_whs2paxinos_tform 
%             for more details.
%
% OUTPUT:
%   ML_whs  - Column vector of coordinates in the medial-lateral direction.
%             Right is positive. Unit in millimeters.
%   AP_whs  - Column vector of coordinates in the anterior-posterior
%             direction. Anterior is positive.  Unit in millimeters.
%   DV_whs  - Column vector of coordinates in the dorsal-ventral direction.
%             Dorsal is positive. Unit in millimeters.

% Get affine transformation:
if nargin==3
    type = 'Anisotropic1';
end
tform = get_whs2paxinos_tform(type); % Augmented matrix form.

% Calculate:
co = tform \ [ML_pax(:) AP_pax(:) DV_pax(:) ones(numel(ML_pax),1)]';
ML_whs = co(1,:)';
AP_whs = co(2,:)';
DV_whs = co(3,:)';

% Reshape:
ML_whs = reshape(ML_whs,size(ML_pax));
AP_whs = reshape(AP_whs,size(AP_pax));
DV_whs = reshape(DV_whs,size(DV_pax));