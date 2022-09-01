classdef Aligner < handle
  properties
    n_ROIs
    scan_data;
    scan_size;
    scan_ROI;
    scan_points;
    scan_pointsC;
    scan_ref;
    scan_resolution;
    atlas_data;
    atlas_size;
    atlas_ROI;
    atlas_points;
    atlas_pointsC;
    atlas_ref;
    atlas_resolution_nominal;
    atlas_resolution_estimated;
    transform_data;
    transform_size;
    transform_ROI;
    transform_points;
    transform_pointsC;
    transform_ref;
    total_transform;
    structure_names;
    segmented_atlas_filename;
    bregma;
    bregmaR;
    bregmaT;
    APvec;
    APvecR;
    APvecT;
    MLvec;
    MLvecR;
    MLvecT;
    DVvec;
    DVvecR;
    DVvecT;
    direct_observation;
  end % properties

  methods
    function obj = Aligner( direct_observation, scan_resolution, n_ROIs, structure_names, segmented_atlas, ... 
                                 bregma, lambda, wxorig, atlas_resolution)
      % scan_resolution and direct_observation are not optional
      % but if left empty the caller will pass -1
      if scan_resolution <= 0.0
        scan_resolution = 0.0343; % millimeters
      end
      obj.scan_resolution    = scan_resolution;
      obj.direct_observation = direct_observation;

      if direct_observation
        % bregma directly observed; skip the whole atlas thing
        disp('Creating Aligner object based on direct landmark observation using the last 4 ROIs!');
        obj.n_ROIs = 4; % default for direct observation

      else % will use the atlas
        disp(['Creating Aligner object with the last ' num2str( n_ROIs) ' ROIs!']);
        obj.n_ROIs = n_ROIs;

        m = load( structure_names);
        obj.structure_names = m.structure_names;
        obj.segmented_atlas_filename = segmented_atlas; % we can load later

        %% check if will use arguments or defaults
        % default values refer to Waxholm Atlas
        if ~exist('bregma', 'var')
          bregma = [246, 653, 440]; % sagittal, coronal, horizontal (adapted for v1.01)
        end

        if ~exist('lambda', 'var')
          lambda = [244, 442, 464]; % sagittal, coronal, horizontal (adapted for v1.01)
        end

        if ~exist('wxorig', 'var')
          wxorig = [244, 623, 248]; % sagittal, coronal, horizontal (adapted for v1.01)
        end

        if ~exist('atlas_resolution', 'var')
          atlas_resolution = 0.039; % millimeters, corrected
        end
        obj.bregma = bregma;
        obj.atlas_resolution_nominal   = atlas_resolution;
        obj.atlas_resolution_estimated = atlas_resolution;

        %% get vectors
        % anterior-posterior
        lb = bregma - lambda;
        obj.APvec = lb / norm( lb);

        % medial-lateral
        lo = wxorig - lambda;
        ml = cross( lo, obj.APvec);
        obj.MLvec = ml / norm( ml);

        % dorsal-ventral
        dv = cross( obj.APvec, -obj.MLvec);
        obj.DVvec = dv / norm( dv); % a bit redundant

        % now that we have all (pixel) unit vectors
        obj.APvec = (obj.APvec/atlas_resolution) + bregma; % 1mm forward
        obj.MLvec = (obj.MLvec/atlas_resolution) + bregma; % 1mm leftward
        obj.DVvec = (obj.DVvec/atlas_resolution) + bregma; % 1mm upward
      end % if direct observation
    end % Aligner
    
    
    function read_scan_data( obj, data, ROI_data)
      %% before this block: prepare and export scan data
      obj.scan_size = size( data);
      assert( isequal( obj.scan_size, size( ROI_data)));

      disp('Reading ROIs! This takes a while...');
      real_scan_points = getRoiCGs( ROI_data);
      % landmarks are in last obj.n_ROIs ROIs
      real_scan_points = real_scan_points( end-obj.n_ROIs+1 : end, :);

      % now's the time to check for direct observation
      if obj.direct_observation
        % calculate vectors and skip the rest
        obj.bregmaT = real_scan_points(1, :);
        lb = obj.bregmaT - real_scan_points(2, :); % bregma - lambda
        obj.APvecT = lb / norm( lb);

        % medial-lateral
        ml = real_scan_points(4, :) - real_scan_points(3, :); % right - left
        obj.MLvecT = ml / norm( ml);

        % double check that the ap-ml angle is very close to perpendicular
        if dot( obj.APvecT, obj.MLvecT) <= 0.0175 % one degree margin %0.0349 % two degree margin
          disp('Landmarks appear consistent.');
        else
          disp('Warning! AP and ML directions as given by the landmarks don''t appear to be perpendicular!');
          disp(['Calculated angle: ' num2str( acosd( dot( obj.APvecT, obj.MLvecT))) ' degrees.']);
        end % angle check

        % dorsal-ventral
        dv = cross( obj.MLvecT, obj.APvecT);
        obj.DVvecT = dv / norm( dv); % a bit redundant

        % now that we have all (pixel) unit vectors
        obj.APvecT = (obj.APvecT/obj.scan_resolution) + obj.bregmaT; % 1mm forward
        obj.MLvecT = (obj.MLvecT/obj.scan_resolution) + obj.bregmaT; % 1mm leftward
        obj.DVvecT = (obj.DVvecT/obj.scan_resolution) + obj.bregmaT; % 1mm upward
      else % I will use the atlas
        obj.scan_data = data;
        obj.scan_ROI = ROI_data;

        obj.scan_points = real_scan_points; % Aligner.ind_coord_swap( real_scan_points);
        % swapping is necessary for the rotations/translations to work out
        % but I will wait until it's actually needed and record real coordinates
        obj.scan_pointsC = [ obj.scan_points; mean( obj.scan_points)]; % copy with centroid
        obj.scan_ref = imref3d( obj.scan_size, 1, 1, 1); % setting pixel extent in world for the math to work out
      end % if direct observation
    end % read_scan_data
    
    
    function read_atlas_data( obj, data, ROI_data)
      %% before this block: prepare and export atlas data 
      obj.atlas_data = data;
      obj.atlas_size = size( data);

      assert( isequal( obj.atlas_size, size( ROI_data)));
      obj.atlas_ROI = ROI_data;
      
      if ~isequal( obj.scan_size, obj.atlas_size)
        %% if resize is necessary, always resize atlas, not scan
        % in order not to merge/lose electrodes
        disp('Atlas size and scan size don''t match! Resizing...');
        
        % check estimates of actual volume covered and
        % update internal resolution estimate
        disp(['Scan volume: '  num2str( (obj.scan_size+1)  * obj.scan_resolution)]);
        disp(['Atlas volume: ' num2str( (obj.atlas_size+1) * obj.atlas_resolution_nominal)]);

        pixel_ratio = obj.scan_size ./ obj.atlas_size;
        obj.atlas_resolution_estimated = obj.atlas_resolution_nominal / mean( pixel_ratio);
        % mean( pixel_ratio) approximately preserves total volume only
        % atlas_aniso_res = obj.atlas_resolution_nominal * pixel_ratio.^-1;
        % preserves everything instead but is anisotropic

        % first scale bregma and vectors accordingly
        obj.bregmaR = obj.bregma .* pixel_ratio;
        obj.MLvecR  = obj.MLvec  .* pixel_ratio;
        obj.APvecR  = obj.APvec  .* pixel_ratio;
        obj.DVvecR  = obj.DVvec  .* pixel_ratio;

        % now resize
        obj.atlas_data = imresize3( obj.atlas_data, obj.scan_size, 'nearest');
        obj.atlas_ROI  = imresize3( obj.atlas_ROI,  obj.scan_size, 'nearest');
        % 'linear' or 'cubic' will destroy ROIs
        
        obj.atlas_size = size( obj.atlas_data); % new size
        assert( isequal( obj.atlas_size, size( obj.atlas_ROI)));

        % if checking is necessary: duplicate and export here
        disp('Don''t forget to re-import from obj.atlas_data and obj.atlas_ROI!');
        %imlook4d_Cdata = atlas_data;
        %imlook4d_ROI = atlas_ROI;
        % then import and check
      end % if size

      % now get points
      disp('Reading ROIs! This takes a while...');
      real_atlas_points = getRoiCGs( obj.atlas_ROI);
      % landmarks are in last obj.n_ROIs ROIs
      real_atlas_points = real_atlas_points( end-obj.n_ROIs+1 : end, :);
      obj.atlas_points = real_atlas_points; % Aligner.ind_coord_swap( real_atlas_points);
      % swapping is necessary for the translations/rotations to work out
      % but I will wait until it's actually needed and record real coordinates
      obj.atlas_pointsC = [ obj.atlas_points; mean( obj.atlas_points)]; % copy with centroid

      obj.atlas_ref = imref3d( obj.atlas_size, 1, 1, 1); % setting pixel extent in world for the math to work out
    end % read_atlas_data
    
    
    function calculate_transform( obj)
      import electroslice.*;
      
      %% first plot points to check
      Aligner.plot_compare_landmarks( obj.atlas_pointsC, obj.scan_pointsC, ...
                                        'Atlas intrinsic', 'Scan intrinsic');

      %% calculating rotation and stretch
      % calculate both point sets centered on origin
      % needed to calculate the transform
      origin_translation_atlas = -obj.atlas_pointsC(end, :);
      origin_translation_scan  = -obj.scan_pointsC(end, :); % last line is centroid

      atlas_points_origin  = obj.atlas_points  + origin_translation_atlas;
      atlas_points_originC = obj.atlas_pointsC + origin_translation_atlas;

      scan_points_origin   = obj.scan_points   + origin_translation_scan;
      scan_points_originC  = obj.scan_pointsC  + origin_translation_scan;
      Aligner.plot_compare_landmarks( atlas_points_originC, scan_points_originC, ...
                                        'Atlas centered', 'Scan centered');

      % calculate svd transform
      [U, S, V, t] = Aligner.svd_transform( atlas_points_origin, scan_points_origin);
      disp( ['Singular Values (degree of asymmetry): ' num2str( [S(1,1) S(2,2) S(3,3)])])
      disp( ['Residual translation should not be too big! norm(t) = ' num2str( norm(t))])

      % composing transformations
      % calculate translations from atlas to origin, then from origin to scan
      aTsvd = [ U*S*V' [0;0;0]; [0 0 0 1]];
      aTatlas_to_origin = [ eye(3) [0;0;0];  origin_translation_atlas 1];
      aTorigin_to_scan  = [ eye(3) [0;0;0]; -origin_translation_scan  1];
      pre_swap = aTatlas_to_origin * aTsvd * aTorigin_to_scan; % transform in canonical coordinates; order is important
      aTrotate = affine3d( Aligner.ind_coord_swap(Aligner.ind_coord_swap( pre_swap)' )' ); % record the swapped transform or it fails
      % the double swap+transpose exchanges the first
      % two rows and the first two columns independently

      %% apply svd transform to the atlas
      disp('Applying transform...');
      [obj.transform_data, obj.transform_ref] = imwarp( obj.atlas_data, ...
                   obj.atlas_ref, aTrotate, 'nearest', 'OutputView', obj.scan_ref);
      [obj.transform_ROI, transform_ROI_ref] = imwarp( obj.atlas_ROI, ...
                   obj.atlas_ref, aTrotate, 'nearest', 'OutputView', obj.scan_ref);

      % for completeness
      obj.transform_size = size( obj.transform_data);

      % check that ROIs still exist:
      disp('Reading ROIs! This takes a while...');
      real_transform_points = getRoiCGs( obj.transform_ROI);
      % landmarks are in last obj.n_ROIs ROIs
      real_transform_points = real_transform_points( end-obj.n_ROIs+1 : end, :);
      disp( "Transformation completed.")
      disp( ['Number of ROIs: '  num2str( numel( real_transform_points(:, 1)))])
      disp( ['Number of NaNs in ROIs: ' num2str( sum( isnan( real_transform_points(:))))])

      % plot results
      obj.transform_points = real_transform_points; % Aligner.ind_coord_swap( real_transform_points);
      % again I'm sticking with unswapped points until it's necessary
      obj.transform_pointsC = [ obj.transform_points; mean( obj.transform_points)];

      Aligner.plot_compare_landmarks( obj.transform_pointsC, obj.scan_pointsC, ...
                                        'Transformed Atlas', 'Scan');
      obj.total_transform = aTrotate;
      disp('Results in obj.transform_data and obj.total_transform.');

      %% clean up and return
      obj.scan_ROI  = [];
      obj.atlas_ROI = [];
      obj.transform_ROI = [];
      disp('Done!');
    end % calculate transform


    function success = refine_alignment( obj, hours)
      import electroslice.*;

      success = false;
      if obj.direct_observation
        disp('This routine is only relevant if you are registering the scan to an Atlas!');
        return;
      elseif isempty( obj.total_transform)
        disp('You need to identify landmarks and calculate the initial alignment first!');
        return;
      end

      disp('Creating Optimizer object to refine Atlas alignment...');
      opt = Optimizer( obj);
      disp(['Running refinement for ' num2str(hours) ' hours...']);
      [resx, resqual] = opt.try_ga( hours*3600);

      if resqual < opt.initial_qual
        disp(['Refinement resulted in ' ...
          num2str(100*(abs(resqual)-abs(opt.initial_qual))/abs(opt.initial_qual)) ...
          '% relative improvement in objective function. Recording new transform...']);
        obj.total_transform = opt.rebuild_transform( opt.aTmap, resx', opt.scale_factor);
      else
        disp('Refinement did not result in improvement of objective function!');
      end % if improvement

      % ready to return!
      success = true;
      disp('Done!');
    end % refine_alignment


    function new_atlas = transform_segmented_atlas( obj, segmented_atlas)
      % transform atlas
      % assumes the reference frame for the segmented atlas
      % is the same as for the actual atlas image
      disp('Transforming the segmented atlas according to obj.total_transform...');
      [new_atlas, new_ref] = imwarp( segmented_atlas, obj.atlas_ref, ...
               obj.total_transform, 'nearest', 'OutputView', obj.scan_ref);
    end % transform_segmented_atlas


    function names = internal_get_names( obj, wires)
      % read segmented atlas
      l = load( obj.segmented_atlas_filename);

      % sizes will almost certainly be different;
      % resize segmented atlas to scan size
      disp('Resizing segmented atlas...');
      l.segmented_atlas = imresize3( l.segmented_atlas, obj.scan_size, 'nearest');

      % transform in helper object
      disp('Transforming segmented atlas...');
      new_atlas = obj.transform_segmented_atlas( l.segmented_atlas);
      clear l;
      % now just read the names
      disp('Reading structure names!');
      names = {};
      for wr = wires
        tip = wr.get_tip();
        value = new_atlas( tip(1), tip(2), tip(3) );
        names{end+1} = obj.structure_names( value);
      end % for wr
      return
    end % get_names
    
    
    function coords = internal_get_stereo_coords( obj, wires)
      import electroslice.*;

      if ~obj.direct_observation
        % transform bregma and unit vectors
        % double swap is needed because the recorded transform is swapped
        obj.bregmaT = Aligner.ind_coord_swap( [obj.bregmaR 1]) * obj.total_transform.T;
        obj.MLvecT  = Aligner.ind_coord_swap(  [obj.MLvecR 1]) * obj.total_transform.T;
        obj.APvecT  = Aligner.ind_coord_swap(  [obj.APvecR 1]) * obj.total_transform.T;
        obj.DVvecT  = Aligner.ind_coord_swap(  [obj.DVvecR 1]) * obj.total_transform.T;

        % swap back
        obj.bregmaT = Aligner.ind_coord_swap( obj.bregmaT(1:3));
        obj.MLvecT  = Aligner.ind_coord_swap(  obj.MLvecT(1:3));
        obj.APvecT  = Aligner.ind_coord_swap(  obj.APvecT(1:3));
        obj.DVvecT  = Aligner.ind_coord_swap(  obj.DVvecT(1:3));
      end % if not direct_observation
      % update resolution estimate
      % after alignment resolution is identical to the scan
      % if using direct observation, naming is just legacy
      obj.atlas_resolution_estimated = obj.scan_resolution;
      
      %% update coords
      coords = {};
      for wr = wires
        vec = wr.get_tip() - obj.bregmaT;
        MLc = (dot( vec, obj.MLvecT - obj.bregmaT) / norm( obj.MLvecT - obj.bregmaT)) * obj.atlas_resolution_estimated;
        APc = (dot( vec, obj.APvecT - obj.bregmaT) / norm( obj.APvecT - obj.bregmaT)) * obj.atlas_resolution_estimated;
        DVc = (dot( vec, obj.DVvecT - obj.bregmaT) / norm( obj.DVvecT - obj.bregmaT)) * obj.atlas_resolution_estimated;
        coords{end+1} = [ MLc, APc, DVc];
      end
      return
    end
  end % methods


  methods (Static)
    % auxiliary function for testing of registration procedures with landmarks
    function plot_compare_landmarks( pointsA, pointsB, legendA, legendB)
      figure;
      scatter3( pointsA(:,1), pointsA(:,2), pointsA(:,3));
      hold on;
      scatter3( pointsB(:,1), pointsB(:,2), pointsB(:,3));
      for i = 1 : size( pointsB',2)
        line( [pointsB(i,1), pointsA(i,1)], ...
              [pointsB(i,2), pointsA(i,2)], ...
              [pointsB(i,3), pointsA(i,3)], 'Color','black');
      end

      % making it extra clear that the labels need to be swapped too
      labels = 'xyz';
      xlabel( labels(1)); ylabel( labels(2)); zlabel( labels(3));
      hold off;
      legend( legendA, legendB);
      clear i;
    end % plot_compare_landmarks


    % auxiliary function for sending subscript indexes and world coordinates
    % basically just swapping x and y
    % could also use worldToSubscript but I don't want to lose track or carry
    % the reference around every time
    function new_points = ind_coord_swap( points)
      % size independent
      swap_mat = eye( size( points, 2));
      swap_mat(1:2, 1:2) = [0 1;
                           1 0];
      new_points = points * swap_mat;
    end % ind_coord_swap


    % expects as inputs Nx3 matrices of 3D points.
    function [U,S,V,t] = svd_transform(A, B)
      assert( nargin == 2)
      assert( isequal(size(A), size(B)) )

      centroid_A = mean(A);
      centroid_B = mean(B);

      % we expect both sets of points to be centered in the origin
      assert( isequal( fix( centroid_A), [0 0 0]) );
      assert( isequal( fix( centroid_B), [0 0 0]) );

      H = (A'*A) \ (A'*B);
      [U,S,V] = svd(H);

      R = V*U';
      if det(R) < 0 % check for reflection
        V(:,3) = V(:,3) * -1; % last eigenvector fixed by hand
        R = V*U';
      end
      t = -R*centroid_A' + centroid_B';
    end % svd_transform
    
    
    function [flip_data, flip_ROI] = helper_orientation( data, ROI)
      flip_data = permute( data, [2 1 3]);
      flip_ROI  = permute( ROI,  [2 1 3]);
    end % helper_orientation
    
    
    function [flip_data, flip_ROI] = helper_flip_pcb( data, ROI)
      flip_data = flip( data, 3);
      flip_ROI  = flip( ROI,  3);
    end % helper_flip_pcb
    
    
    function [flip_data, flip_ROI] = helper_flip_snout( data, ROI)
      flip_data = flip( data, 2);
      flip_ROI  = flip( ROI,  2);
    end % helper_flip_snout
    
    
    function [flip_data, flip_ROI] = helper_flip_lateral( data, ROI)
      flip_data = flip( data, 1);
      flip_ROI  = flip( ROI,  1);
    end % helper_flip_lateral


    function [flip_data, flip_ROI] = helper_cycle( data, ROI)
      flip_data = permute( data, [3 1 2]);
      flip_ROI  = permute( ROI,  [3 1 2]);
    end % helper_cycle
  end % methods (Static)
end % class
