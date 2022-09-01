function [tform_scale, tform_rot, tform_trans] = decompose_affine_tform(tform)
%DECOMPOSE_AFFINE_TFORM   Decompose affine transformation.
%
% Decompose affine transformation into scaling, rotation and translation.
% The solution assumes
%   1. that the translation is done before the scaling but after the 
%      rotation, like
%      t = t_scale * t_trans * t_rot
%   2. that the scaling is uniform.
%   3. no shearing.
%
% SYNTAX:
%   [tform_scale, tform_rot, tform_trans] = decompose_affine_tform(tform)
%
% INPUT:
%   tform   - Affine transformation on augmented matrix form.
%
% OUTPUT:
%   tform_scale - Scale transformation.
%   tform_rot   - Rotation transformation.
%   tform_trans - Translation transformation.

% Check input:
if size(tform,1)~=4 || size(tform,2)~=4
    error('Input matrix must be 4-by-4.')
end
if ~all(tform(4,1:3)==0)
    error('Shearing is not allowed.')
end
dummy = tform(1:3,1:3)*tform(1:3,1:3)';
dummy = dummy.^2;
if sum(dummy(:))-trace(dummy) > 1e-5
    error('This is probably not a uniform scaling.')
end

% Calculate scaling matrix:
s = det(tform)^(1/3);            % The determinant equals the volume scaling.
tform_scale = s*eye(4); tform_scale(4,4) = 1;
tform_rigid = tform_scale\tform; % Assuming t = t_scale * t_rigid.

% Split t_rigid into t_rot and t_trans:
tform_rot = tform_rigid; tform_rot(:,4) = [0;0;0;1];
tform_trans = eye(4); tform_trans(1:3,4) = tform_rigid(1:3,4); % Assuming tform_rigid=tform_trans*tform_rot