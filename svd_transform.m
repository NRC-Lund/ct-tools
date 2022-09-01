function [U,S,V,t,T] = svd_transform(A, B)
% SYNTAX:
%   [U,S,V,t,T] = svd_transform(A, B)


% expects as inputs Nx3 matrices of 3D points.
assert( nargin == 2)
assert( isequal(size(A), size(B)) )

centroid_A = mean(A);
centroid_B = mean(B);

% % we expect both sets of points to be centered in the origin
% assert( isequal( fix( centroid_A), [0 0 0]) );
% assert( isequal( fix( centroid_B), [0 0 0]) );

A = A - centroid_A;
B = B - centroid_B;

H = (A'*A) \ (A'*B); % same as inv(A'*A) * (A'*B)
[U,S,V] = svd(H);

R = V*U';
if det(R) < 0 % check for reflection
    V(:,3) = V(:,3) * -1; % last eigenvector fixed by hand
    %R = V*U';
end
t = (-centroid_A*(U*S*V') + centroid_B)';

% Affine transformation, represented as an augmented matrix:
T = [U*S*V' t; [0 0 0 1]];
