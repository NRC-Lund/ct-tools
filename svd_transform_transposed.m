function [U,S,V,t,T] = svd_transform_transposed(A, B)
% SYNTAX:
%   [U,S,V,t,T] = svd_transform_transposed(A, B)


% expects as inputs 3xN matrices of 3D points.
assert( nargin == 2)
assert( isequal(size(A), size(B)) )

centroid_A = mean(A, 2);
centroid_B = mean(B, 2);

% % we expect both sets of points to be centered in the origin
% assert( isequal( fix( centroid_A), [0 0 0]) );
% assert( isequal( fix( centroid_B), [0 0 0]) );

A = A - centroid_A;
B = B - centroid_B;

H = (B*A') / (A*A'); % same as (B*A') * inv(A*A')
[U,S,V] = svd(H);

R = V*U';
if det(R) < 0 % check for reflection
    V(:,3) = V(:,3) * -1; % last eigenvector fixed by hand
    %R = V*U';
end
t = -(U*S*V')*centroid_A + centroid_B;

% Affine transformation, represented as an augmented matrix:
T = [U*S*V' t; [0 0 0 1]];
