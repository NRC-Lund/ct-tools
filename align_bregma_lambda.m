function [tform_scale, tform_roll, tform_align, tform_trans] = ...
    align_bregma_lambda(Bregma, Lambda, Left, Right, Midline, LambdaBregmaDistance)
% SYNTAX:
%   [tform_scale, tform_roll, tform_align, tform_trans] = ...
%       align_bregma_lambda(Bregma, Lambda, Left, Right, Midline, LambdaBregmaDistance)
%
% INPUT:
%   Bregma  - [ML,AP,DV]
%   Lambda  - [ML,AP,DV]
%   Left    - [ML1,AP1,DV1; ...]
%   Right   - [ML1,AP1,DV1; ...]
%   Midline - [ML1,AP1,DV1; ...]
%   LambdaBregmaDistance - scalar
%   
%
% OUTPUT:
%   tform   - Affine transformation in augmented matrix form.

% Translation to bregma:
co = [Bregma 1; Lambda 1]';
tform_trans = [1 0 0 -Bregma(1); ...
    0 1 0 -Bregma(2); ...
    0 0 1 -Bregma(3); ...
    0 0 0 1];
co = tform_trans * co;

% Rotation to point Lambda in negative AP direction:
A = co(1:3,2);
B = [0 -1 0]';
A = A / norm(A,2);
B = B / norm(B,2);
DotAB = dot(A,B);
CrossAB = cross(A,B);
G = [DotAB -norm(CrossAB) 0;...
      norm(CrossAB) DotAB  0;...
      0	0 1];
F = [A (B-DotAB*A)/norm(B-DotAB*A) -CrossAB];
tform_align = F*G/F;
tform_align = [tform_align [0;0;0]; 0 0 0 1];
co = tform_align * co;

% Roll to allign ML-AP plane:
if nargin==2
    tform_roll = eye(4);
else
    % Align:
    tform = affine3d((tform_align*tform_trans)');
    [X,Y,Z] = tform.transformPointsForward(Left(:,1),Left(:,2),Left(:,3));
    B_left = [X Y Z];
    [X,Y,Z] = tform.transformPointsForward(Right(:,1),Right(:,2),Right(:,3));
    B_right = [X Y Z];
    [X,Y,Z] = tform.transformPointsForward(Midline(:,1),Midline(:,2),Midline(:,3));
    B_mid = [X Y Z];
    % Find roll:
    f = @(a) norm(B_mid(:,1)*cos(a) + B_mid(:,3)*sin(a),2) + ... % Deviation from midline
        norm((B_left(:,1)-B_right(:,1))*sin(a)-(B_left(:,3)-B_right(:,3))*cos(a),2); % Left-right difference in DV
    a = fminbnd(f,-pi/2,pi/2); % Angle that minimizes f.
    % Rotate:
    tform_roll = [cos(a) 0 sin(a) 0; 0 1 0 0; -sin(a) 0 cos(a) 0; 0 0 0 1];
end

% Scale:
S = LambdaBregmaDistance / norm(Lambda-Bregma,2);
tform_scale = [S 0 0 0; ...
    0 S 0 0; ...
    0 0 S 0; ...
    0 0 0 1];
co = tform_scale * co;

% Checks:
tol = 1e-10;
if co(2,2)>0
    error('Lambda AP should be negative.')
end
if co(1,2)>tol
    error('Lambda is not on the midline.')
end
if co(3,2)>tol
    error('Lambda is not on DV=0.')
end
