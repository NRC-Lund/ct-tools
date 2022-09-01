function tform = get_whs2paxinos_tform(type)
% Returns the affine transformation needed to convert coordinates from the 
% Waxholm atlas (WHS SD rat v1.01) to the Paxinos atlas (Paxinos & Watson,
% 2007). The coordinate system is defined as Right, Anterior, Superior with
% units in mm.
%
% SYNTAX:
%   tform = get_whs2paxinos_tform(type)
%
% INPUT:
%   type    - String specifying the transformation type. Possible values
%             are 'Papp', 'Modified Papp', 'Anisotropic1' and 'Legacy'.
%
% OUTPUT:
%   tform   - Affine transformation on augmented matrix form.
%
% NOTES ON TRANSFORMATION TYPES:
%   'Papp'          - Following Papp et al., 2016 (https://doi.org/10.3389/fninf.2016.00011).
%                     I.e. a translation to bregma, followed by a uniform
%                     scaling of 1.057 and a -4.085 degrees rotation around
%                     the ML axis.
%   'Modified Papp' - Like 'Papp', but with modified rotations to bring
%                     bregma and lambda to DV=0.0000 mm, and bregma, lambda 
%                     and WHS_origo to ML=0.0000 mm. The rotations around
%                     each axis are ML=-6.489, DV=-0.5395 and AP=0.5675 (in
%                     degrees) and are done in that order.
%   'Anisotropic1'  - Like 'Modified Papp', but with anisotropic scaling 
%                     ML=0.9412, AP=1.0491, DV=0.9244.
%   'Legacy'        - Legacy transformation calculated using Horn's
%                     quaternion-based method.
%
%   Table comparing the scaling in the AP, DV and ML directions:
%
%                     bregma-lambda | bregma-WHS_origo | left_dcw-right_dcw
%   Paxinos             8.7033              7.0045              12.8
%   Papp                8.7686              8.0241              14.3752
%   Modified Papp       8.7033              7.9644              14.2682
%   Anisotropic1        8.7033              7.0196              12.8001
%   Legacy              8.3892              7.6770              13.7522
%
%   (Ideally, the distances should be as close as possible to the Paxinos
%   values.)


debugFlag = 0;
switch type
    case 'Legacy'
        if debugFlag
            % WHS coordinates (ML, AP, DV): (Obtained from https://www.nitrc.org/docman/view.php/1081/2095/Coordinates_v1-v1.01.pdf)
            L1 = [0.0781250 1.1718750 7.5; .... % bregma
                0 -7.0703125 8.4375; ...        % lambda
                0 0 0];                         % WHS origo
            % Paxinos coordinates (ML, AP, DV):
            L2 = [0 0 0; ...                    % bregma
                0 -9 0; ...                     % lambda
                0 -0.3 -7];                     % WHS origo
            % Calculate affine transformation:
            [regParams,Bfit,ErrorStats] = absor(L1', L2', ...
                'doTrans', true, 'doScale', true, 'weights', [0.6 0.2 0.2]);
            tform = regParams.M; % Augmented matrix form.
        else % Precalculated:
            tform = [1.0111   -0.0106   -0.0089   -0.0000; ...
                     0.0093    1.0023   -0.1341   -0.3183; ...
                     0.0102    0.1340    1.0023   -7.5068; ...
                     0         0         0         1.0000];
        end
    case 'Papp' % Papp et al., 2016 ("Brain-Wide Mapping of Axonal Connections: Workflow for Automated Detection and Spatial Analysis of Labeling in Microscopic Sections")
        if debugFlag
            tform_trans = [1 0 0 -0.0781250; ...
                0 1 0 -1.1718750; ...
                0 0 1 -7.5; ...
                0 0 0 1]; % Translation to bregma as defined in WHS v1.01.
            s = 1.057;    % Uniform scaling according to Papp.
            tform_scale = [s 0 0 0; ...
                0 s 0 0; ...
                0 0 s 0; ...
                0 0 0 1];
            a = -4.085/180*pi; % Pitch rotation according to Papp.
            tform_rot = [1 0 0 0; ...
                0 cos(a) sin(a) 0; ...
                0 -sin(a) cos(a) 0; ...
                0 0 0 1];
            tform = tform_rot * tform_scale * tform_trans;
        else % Precalculated:
            tform = [1.057000000000000  0                   0                  -0.082578125000000; ...
                     0                  1.054314656340088  -0.075296782311605  -0.670799120561507; ...
                     0                  0.075296782311605   1.054314656340088  -7.995598339322074; ...
                     0                  0                   0                   1];
        end
    case 'Modified Papp' % Improved rotations to put lambda, bregma and WHS origo at ML=0.
        if debugFlag
            % WHS coordinates (ML, AP, DV):
            L1 = [0.0781250 1.1718750 7.5; ...  % bregma
                0 -7.0703125 8.4375; ...        % lambda
                0 0 0];                         % WHS origo
            tform_trans = [1 0 0 -L1(1,1); ...
                0 1 0 -L1(1,2); ...
                0 0 1 -L1(1,3); ...
                0 0 0 1]; % Translation to bregma.
            s = (9.1-0.3)*9/9.1/norm(L1(1,:)-L1(2,:));    % Scaling to Paxinos bregma-lambda distance.
            tform_scale = [s 0 0 0; ...
                0 s 0 0; ...
                0 0 s 0; ...
                0 0 0 1];
            a_ML = -6.489/180*pi; % Pitch
            tform_pitch = [1 0 0 0; ...
                0 cos(a_ML) sin(a_ML) 0; ...
                0 -sin(a_ML) cos(a_ML) 0; ...
                0 0 0 1];
            a_DV = -0.5395/180*pi; % Yaw
            tform_yaw = [cos(a_DV) sin(a_DV) 0 0; ...
                -sin(a_DV) cos(a_DV) 0 0; ...
                0 0 1 0; ...
                0 0 0 1];
            a_AP = 0.5675/180*pi; % Roll
            tform_roll = [cos(a_AP) 0 -sin(a_AP) 0; ...
                0 1 0 0; ...
                sin(a_AP) 0 cos(a_AP) 0; ...
                0 0 0 1];
            tform = tform_roll * tform_yaw * tform_pitch * tform_scale * tform_trans;
        else % Precalculated:
            tform = [1.049035324444181  -0.010989117314430  -0.009208311479476  -0.000015676773282; ...
                     0.009878546942000   1.042365883617422  -0.118559877750217  -0.333095198217384; ...
                     0.010390765610293   0.118462101967175   1.042372019971235  -7.957424704090347; ...
                     0                   0                   0                  1];
        end
    case 'Anisotropic1'
        if debugFlag
            % WHS coordinates (ML, AP, DV):
            L1 = [0.0781250 1.1718750 7.5; ...  % bregma
                  0 -7.0703125 8.4375; ...      % lambda
                  0 0 0];                       % WHS origo
            % Paxinos coordinates (ML, AP, DV):
            L2 = [0 0 0; ...                    % bregma
                  0 -(9.1-0.3)*9/9.1 0; ...     % lambda (The bregma-lambda distance is 9.1-0.3 according to Table 1 in Paxinos, but in the figures the bregma-interaural distance is rescaled from 9.1 to 9)
                  0 -0.5 -7];                   % WHS origo
            % Distances:
            Lambda_Bregma_Wax = norm(L1(1,:)-L1(2,:));
            Lambda_Bregma_Pax = norm(L2(1,:)-L2(2,:));
            AC_Bregma_Wax = norm(L1(1,:)-L1(3,:));
            AC_Bregma_Pax = norm(L2(1,:)-L2(3,:));
            DCW_Wax = 6.8; % The distance between the most lateral edges of the deep cerebral white matter (Wax 6.8,-5.7,1.45).
            DCW_Pax = 6.4; % The distance between the most lateral edges of the deep cerebral white matter (Pax 6.4,-6.48,-5.5).
            % Transformations:
            tform_trans = [1 0 0 -0.0781250; ...
                0 1 0 -1.1718750; ...
                0 0 1 -7.5; ...
                0 0 0 1]; % Translation to bregma.
            s = [DCW_Pax/DCW_Wax ...                    % ML scaling to DCW distances
                Lambda_Bregma_Pax/Lambda_Bregma_Wax ... % AP scaling to Paxinos bregma-lambda distances
                AC_Bregma_Pax/AC_Bregma_Wax];           % DV scaling to AC-bregma distances
            tform_scale = [s(1) 0 0 0; ...
                0 s(2) 0 0; ...
                0 0 s(3) 0; ...
                0 0 0 1];
            a_ML = -6.489/180*pi; % Pitch
            tform_pitch = [1 0 0 0; ...
                0 cos(a_ML) sin(a_ML) 0; ...
                0 -sin(a_ML) cos(a_ML) 0; ...
                0 0 0 1];
            a_DV = -0.5395/180*pi; % Yaw
            tform_yaw = [cos(a_DV) sin(a_DV) 0 0; ...
                -sin(a_DV) cos(a_DV) 0 0; ...
                0 0 1 0; ...
                0 0 0 1];
            a_AP = 0.5675/180*pi; % Roll
            tform_roll = [cos(a_AP) 0 -sin(a_AP) 0; ...
                0 1 0 0; ...
                sin(a_AP) 0 cos(a_AP) 0; ...
                0 0 0 1];
            tform = tform_scale * tform_roll * tform_yaw * tform_pitch * tform_trans;
        else % Precalculated:
            tform = [0.941088583454037  -0.009858326603374  -0.008260767396766  -0.000014063618276; ...
                     0.009878546942000   1.042365883617422  -0.118559877750217  -0.333095198217384; ...
                     0.009155835230371   0.104383019243590   0.918487320525401  -7.011694054243962; ...
                     0                   0                   0                   1];
        end
end