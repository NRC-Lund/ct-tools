classdef Optimizer
  properties
    atlas
    scan
    transform_ref
    scale_factor
    aT
    aTmap
    skull_mask
    initial_qual
  end % properties


  methods
    function obj = Optimizer( a_obj, scale_factor, gaussian_px)
      import electroslice.*;
      
      % default values
      if ~exist('scale_factor', 'var')
        scale_factor = 0.5;
      end
      
      if ~exist('gaussian_px', 'var')
        gaussian_px = fix( scale_factor * 50);
      end
      
      % prepare data (right now just downsample)
      obj.scale_factor = scale_factor;
      obj.atlas = imresize3( a_obj.atlas_data,     scale_factor);
      obj.scan  = imresize3( a_obj.scan_data,      scale_factor);
      t_data    = imresize3( a_obj.transform_data, scale_factor);

      % prepare transforms
      assert( all(size( obj.scan) == round((a_obj.scan_size)*scale_factor)));
      obj.transform_ref = imref3d( round((a_obj.scan_size)*scale_factor), 1, 1, 1); % resized ref
      obj.aT    = a_obj.total_transform;
      obj.aTmap = Optimizer.decompose_transform( obj.aT);
      obj.aTmap(22:24) = obj.aTmap(22:24) * scale_factor; % rescale translation

      % smooth and threshold atlas to create skull mask
      disp('Calculating filter to isolate skull from Atlas...');
      obj.skull_mask  = Optimizer.gaussian3d( obj.atlas, gaussian_px);
      skull_tdata     = Optimizer.gaussian3d( t_data,    gaussian_px);

      % note:
      skull_tdata2 = imwarp( obj.skull_mask, obj.transform_ref, ...
        obj.rebuild_transform(obj.aTmap), 'nearest', 'OutputView', obj.transform_ref);
      skull_tdata2( skull_tdata2 < 1500) = 0;
      skull_tdata2( skull_tdata2 >= 1500) = 1;
      skull_tdata2 = logical( skull_tdata2);
      % results in a very similar mask (~1% of pixels are different)
      % and superposition index is ~0.75, and saves ~20min of computation.
      % Maybe use that? difference in quality measure is ~2%

      obj.skull_mask( obj.skull_mask < 2000)  = 0;
      obj.skull_mask( obj.skull_mask >= 2000) = 1;
      obj.skull_mask = logical( obj.skull_mask);
      skull_tdata( skull_tdata < 2000) = 0;
      skull_tdata( skull_tdata >= 2000) = 1;
      skull_tdata = logical( skull_tdata);
      
      obj.atlas( obj.skull_mask) = ( max( obj.atlas(:)) - obj.atlas( obj.skull_mask)).^4;
      obj.atlas = obj.atlas .* 10^-14; % keep intensities under control

      t_data( skull_tdata) = ( max( t_data(:)) - t_data( skull_tdata)).^4;
      t_data = t_data .* 10^-14;
      
      % threshold scan to remove screws and metal
      num_metal = numel(find( obj.scan > 7000));
      noise = rand( num_metal, 1) .* 500;
      obj.scan( obj.scan > 7000) = noise;

      % makes more sense to keep the transformed mask if i'm not warping
      obj.skull_mask = skull_tdata;
      
      % get initial quality value
      obj.initial_qual = Optimizer.superposition_index( t_data, obj.scan, skull_tdata);
    end % Optimizer


    function [res1, res2] = try_fmincon( obj, N, initial_guess)
      import electroslice.*;

      % set bounds
      % reminder: axis perturb direction 1, axis perturb angle 1, delta angle 1, 
      %           axis perturb direction 2, axis perturb angle 2, delta angle 2, 
      %           delta scale 1 2 3, delta translation 1 2 3.
      lb = [ -eps;-eps;-0.15;  -eps;-eps;-0.15;  -0.1;-0.1;-0.1;  -8;-8;-8];
      ub = [ 2*pi;0.15; 0.15;  2*pi;0.15; 0.15;   0.1; 0.1; 0.1;   8; 8; 8];
      
      % initial guess is no perturbation unless supplied
      if ~exist('initial_guess', 'var')
        initial_guess = zeros(12,1);
      end

      % go
      [res1, res2] = fmincon( @(nn) obj.superposition_index( ...
                       imwarp( obj.atlas, obj.transform_ref, ...
                       Optimizer.rebuild_transform( obj.aTmap, nn),...
                       'nearest', 'OutputView', obj.transform_ref), obj.scan,...
                       obj.skull_mask), ...
                       initial_guess, [],[],[],[], lb, ub, [], ...
                       optimoptions(@fmincon,'MaxIterations',N,...
                       'Display', 'iter'));
                       % obj.skull_mask used to be:
                       %imwarp( obj.skull_mask, obj.transform_ref, ...
                       %Optimizer.rebuild_transform( obj.aTmap, nn), ...
                       %'nearest', 'OutputView', obj.transform_ref), ...
    end % try_fmincon


    function [res1, res2] = try_ga( obj, N)
      import electroslice.*;

      % set bounds
      % reminder: axis perturb direction 1, axis perturb angle 1, delta angle 1, 
      %           axis perturb direction 2, axis perturb angle 2, delta angle 2, 
      %           delta scale 1 2 3, delta translation 1 2 3.
      lb = [ -eps;-eps;-0.15;  -eps;-eps;-0.15;  -0.1;-0.1;-0.1;  -8;-8;-8];
      ub = [ 2*pi;0.15; 0.15;  2*pi;0.15; 0.15;   0.1; 0.1; 0.1;   8; 8; 8];

      % go
      [res1, res2] = ga( @(nn) obj.superposition_index( ...
                       imwarp( obj.atlas, obj.transform_ref, ...
                       Optimizer.rebuild_transform( obj.aTmap, nn' ),...
                       'nearest', 'OutputView', obj.transform_ref), ...
                       obj.scan, ...
                       obj.skull_mask), ...
                       12, [],[],[],[], lb, ub, [], ...
                       optimoptions(@ga,'MaxTime',N, 'PlotFcn', @gaplotscores));
                       % obj.skull_mask used to be:
                       %imwarp( obj.skull_mask, obj.transform_ref, ...
                       %Optimizer.rebuild_transform( obj.aTmap, nn), ...
                       %'nearest', 'OutputView', obj.transform_ref), ...
    end % try_ga
  end % methods


  methods (Static)
    function value = superposition_index( atlas, scan, mask)
      % mutual information; from spm
      % joint histogram
      jh = histcounts2( atlas( mask), scan( mask), 'Normalization', 'probability');

      jh = jh + eps; % avoiding NaN/Inf
      jh = jh / sum( jh(:)); % normalize again
      s1 = sum(jh, 1);
      s2 = sum(jh, 2);

      % mutual information: between zero and one, higher is better!
      jh = jh .* log2( jh ./ (s2*s1));
      value = -1 * sum(jh(:)); % but we return -value to agree with the minimizer
    end % superposition_index


    function aT = rebuild_transform( aTmap, perturb, scale_factor)
      import electroslice.*;

      % rebuild matrices and vector
      if exist('perturb', 'var')
       aTmap([7:9 16:24]) = aTmap([7:9 16:24]) + perturb;
      end % if exist
      
      if exist('scale_factor', 'var')
        aTmap(22:24) = aTmap(22:24) ./ scale_factor;
      end
      
      % apply perturbation
      vec = Optimizer.perturb_in_cone( ...
              aTmap(1:3), aTmap(4:6), aTmap(7), aTmap(8));
      % create rotation matrix
      pU = Optimizer.mat_from_vec_and_angle( vec, aTmap( 9) );
      
      % apply perturbation
      vec = Optimizer.perturb_in_cone( ...
              aTmap(10:12), aTmap(13:15), aTmap(16), aTmap(17));
      % create rotation matrix
      pV = Optimizer.mat_from_vec_and_angle( vec, aTmap( 18) );

      % create scaling and translation
      pS = diag( aTmap( 19:21));
      pD = aTmap( 22:24);

      % compose
      aT = affine3d( [ pU*pS*pV [0;0;0]; pD' 1] ); % notice the lack of transposition here!
    end


    function aTmap = decompose_transform( aT)
      import electroslice.*;

      % get rotation matrix and translation vector
      origT = aT.T;
      origR = origT(1:3, 1:3);
      origD = origT(end, 1:3)';

      % svd decomposition
      [u,s,v] = svd( origR);

      % check for total reflection
      % shouldn't really ever happen
      R = v*u';
      if det(R) < 0
        v(:,3) = v(:,3) * -1; % last eigenvector fixed by hand
      end

      % check for pair of reflections
      if (det(u) < 0) && (det(v) < 0)
        u = u*v;
        v = v*v;
      end
      
      % initialize and fill
      aTmap = zeros(24, 1);
      [aTmap(1:3), aTmap(9)]     = Optimizer.vec_and_angle_from_mat( u);
      aTmap(4:6)   = Optimizer.orthogonal( aTmap(1:3) );   % orthogonal direction
      aTmap(7)   = 0.0; % angle for perturbation direction
      aTmap(8)   = 0.0; % cone angle
      [aTmap(10:12), aTmap(18)]  = Optimizer.vec_and_angle_from_mat( v'); % notice the transposition
      aTmap(13:15) = Optimizer.orthogonal( aTmap(10:12) ); % orthogonal direction
      aTmap(16)  = 0.0; % angle for perturbation direction
      aTmap(17)  = 0.0; % cone angle
      aTmap(19:21)               = diag( s);
      aTmap(22:24)               = origD;
      assert(norm(origR-u*s*v')<0.001)
    end


    function ov = orthogonal( vector)
      vector = vector ./ norm( vector);
      if ismembertol( vector(3), 0)
        vector(3) = 10*eps; % avoid division by zero
      end
      vector = vector ./ norm( vector);

      ov = [1;1; -(vector(1)+vector(2))/vector(3)];
      ov = ov ./ norm( ov);
 
      assert( ismembertol( norm( cross(ov, vector)), 1));
    end


    function pv = perturb_in_cone( vector, orthov, angle_dir, angle_cone)
      import electroslice.*;

      % assumes all vectors are unit vectors
      % angle_dir should be 0 <= angle_dir <= 2pi
      % angle_cone should be small
      
      % rotate perpendicular direction around vector
      rot_orthov = Optimizer.mat_from_vec_and_angle( vector, angle_dir) * orthov;
      rot_orthov = rot_orthov ./ norm( rot_orthov); % why not
      
      % rotate vector around rotated perpendicular direction
      pv = Optimizer.mat_from_vec_and_angle( rot_orthov, angle_cone) * vector;
      pv = pv ./ norm( pv); % why not
    end


    function check = is_rotation_matrix( M)
      check = ( det(M) > 0); % no reflection
      check = ( abs( det(M) - 1.0)  < 0.001) && check;
      check = ( norm(M*M' - eye(3)) < 0.001) && check;
      check = ( norm(M'*M - eye(3)) < 0.001) && check;
    end


    function M = mat_from_vec_and_angle( vec, angle)
      import electroslice.*;

      % imposing all angles between 0 and 2pi
      if angle < 0
        angle = 2*pi + angle;
      end
      % vector defined up to direction sign
      % imposing vector direction is towards positive octant
      if sum( (vec < 0.0)) > 1 % has 2 or 3 negative components
        vec = -vec;
        angle = 2*pi - angle; % change rotation direction to match
      end
      M = cos( angle)*eye(3) + ...
          sin( angle)*[cross( vec,[1 0 0])' cross( vec,[0 1 0])' cross( vec,[0 0 1])'] + ...
          (1-cos( angle))*(vec*vec');
      % make sure it's a rotation matrix
      assert( Optimizer.is_rotation_matrix( M));
    end % mat_from_vec_and_angle


    function [vec, angle] = vec_and_angle_from_mat( M)
      import electroslice.*;

      % make sure it's a rotation matrix
      assert( Optimizer.is_rotation_matrix( M));

      angle = acos( (trace( M)-1)/2 ); % between 0 and pi
      possible_solutions = [angle, 2*pi-angle]; % solutions between 0 and 2pi

      [dummy1, dummy2] = eig( M);

      if all( imag(diag(dummy2)) == 0)
        % if input is diagonal matrix it can happen that all eigenvalues
        % are real; in this case it's better to look for the actual
        % eigenvector whose eigenvalue is one
        if sum( ismembertol(diag(dummy2), 1)) == 1 % only one eigenvalue == 1
          vec = dummy1(:, ismembertol( diag(dummy2), 1));
        else % all eigenvalues==1, it's probably the identity matrix?
          vec = dummy1(:,1); % convention
        end
      else
        vec = dummy1(:, imag(diag(dummy2)) == 0 ); % eigenvalue == 1, safer to detect zero imaginary part
      end

      % vector defined up to direction sign
      % imposing vector direction is towards positive octant
      if sum( (vec < 0.0)) > 1 % has 2 or 3 negative components
        vec = -vec;
      end

      if norm( M-M') > 0.01
        % not symmetrical, I can use this method
        u1 = [ M(3,2)-M(2,3) M(1,3)-M(3,1) M(2,1)-M(1,2)];
        assert( norm(cross( vec, u1)) < 0.01); % parallel
        
        direc = dot(u1, vec)/(norm(u1)*norm(vec)); % norm(vec) is expected to be 1 but you never know
        angle1 = asin( direc*norm( u1)/2); % between -pi/2 and pi/2

        possible_solutions2 = [angle1, pi-angle1];
        if possible_solutions2(1) < 0.0
          possible_solutions2(1) = 2*pi + possible_solutions2(1);
        end

        assert( all( possible_solutions2 >= 0))

        % pick the one that is in both
        angle = possible_solutions( ismembertol( possible_solutions, possible_solutions2));
        if numel(angle) == 0
          % some numerical problems when the angle is pi/2, 3pi/2
          angle = possible_solutions( ismembertol( possible_solutions, possible_solutions2, 1e-8));
          assert( (abs(angle -pi/2) < 0.001) || abs(angle -3*pi/2) < 0.001);
        end
        assert( numel(angle) == 1); % exactly one
      else
        % we need to pick the right solution somehow
        % pick the one that results in the same matrix
        sol = arrayfun( @(nn) ...
          norm(M - Optimizer.mat_from_vec_and_angle( vec, nn)), possible_solutions);
        [~, ind] = min(sol);
        angle = possible_solutions( ind);
      end
    end % vec_and_angle_from_mat


    function testing_plot( argument, type)
      obj = argument.align_helper;
      if type == 1 % scan coordinate system
        figure()
        title('Wire tips in Scan coordinates with transformed direction vectors');
        hold on;
        scatter3( obj.bregmaT(2), obj.bregmaT(1), obj.bregmaT(3));
        text( obj.MLvecT(2), obj.MLvecT(1), obj.MLvecT(3), 'ML');
        text( obj.DVvecT(2), obj.DVvecT(1), obj.DVvecT(3), 'DV'); 
        text( obj.APvecT(2), obj.APvecT(1), obj.APvecT(3), 'AP');
        for wr = argument.wires
          coords = wr.get_tip();
          scatter3( coords(2), coords(1), coords(3));
        end
        grid();
        quiver3( obj.bregmaT(2), obj.bregmaT(1), obj.bregmaT(3), obj.MLvecT(2)-obj.bregmaT(2), obj.MLvecT(1)-obj.bregmaT(1), obj.MLvecT(3)-obj.bregmaT(3));
        quiver3( obj.bregmaT(2), obj.bregmaT(1), obj.bregmaT(3), obj.DVvecT(2)-obj.bregmaT(2), obj.DVvecT(1)-obj.bregmaT(1), obj.DVvecT(3)-obj.bregmaT(3));
        quiver3( obj.bregmaT(2), obj.bregmaT(1), obj.bregmaT(3), obj.APvecT(2)-obj.bregmaT(2), obj.APvecT(1)-obj.bregmaT(1), obj.APvecT(3)-obj.bregmaT(3));
        hold off;
      else % stereo coordinate system
        figure()
        title('Wire tips in Stereo coordinates');
        hold on;
        scatter3( 0, 0, 0);
        text( 1, 0, 0, 'ML'); text( 0, 0, 1, 'DV'); text( 0, 1, 0, 'AP');
        for wr = argument.wires
          coords = wr.stereo_coords;
          scatter3( coords(2), coords(1), coords(3));
       end
       grid();
       quiver3( 0,0,0, 1, 0, 0);
       quiver3( 0,0,0, 0, 0, 1);
       quiver3( 0,0,0, 0, 1, 0);
       hold off;
      end
    end % testing_plot


    function matrix3d = gaussian2d( matrix3d, fwhm)
      fsize=( round(3*fwhm) +1 )/2; % This will be an integer or an integer+0.5
      sigma = fwhm/2.35;

      % Make Gaussian Kernel
      for i = -(fsize-1):(fsize-1)
        for j = -(fsize-1):(fsize-1)
          k = i+fsize;
          l = j+fsize;
          filterHandle(k,l) = exp(-(i^2+j^2)/(2*sigma^2));
        end
      end

      filterHandle=filterHandle/sum(filterHandle(:));  % Normalize to sum 1

      % Loop  slices
      for i=1:size(matrix3d,3)
           temp=conv2(matrix3d(:,:,i),filterHandle,'same');
           matrix3d(:,:,i)=temp;
      end
    end % gaussian2d


    function mat = gaussian3d( mat_in, pixels)
      % FWHM-in-pixels / 2.35
      sigmaX = pixels / 2.35;
      sigmaY = sigmaX;
      sigmaZ = sigmaX;

      c = 5; % Number of sigmas to use in matrix (convolution speed depends on this)
      fsizeX = ceil( c*sigmaX);
      fsizeY = ceil( c*sigmaY);
      fsizeZ = ceil( c*sigmaZ);

      % assuming
      dX = 1; dY = 1; dZ = 1;

      % suppress prints for now
      %disp(['Filtersize sigma  (x,y,z)=(' num2str(sigmaX) ', ' num2str(sigmaY) ', ' num2str(sigmaZ) ') pixels']);
      %disp(['Filtersize pixels (x,y,z)=(' num2str(fsizeX) ', ' num2str(fsizeY) ', ' num2str(fsizeZ) ') pixels']);

      % Gaussian filter matrix
      for i = -(fsizeX-1):(fsizeX-1)  % Loop in steps of 1 pixel
        for j = -(fsizeY-1):(fsizeY-1)
          for k= -(fsizeZ-1):(fsizeZ-1)
            x = i*dX; % position in mm
            y = j*dY; % position in mm
            z = k*dZ; % position in mm
            % Calculate formula using positions in mm, place in ix=1,2,.., iy=1,2,..., iz=1,2,...
            %filterKernel(i+fsizeX, j+fsizeY, k+fsizeZ) = exp( -(x^2+y^2+z^2)/(2*sigmaX^2+2*sigmaY^2+2*sigmaZ^2)  );
            
            filterKernel(i+fsizeX, j+fsizeY, k+fsizeZ) = ...
              exp( -(x^2)/(2*sigmaX^2) ) * ...
              exp( -(y^2)/(2*sigmaY^2) ) * ...
              exp( -(z^2)/(2*sigmaZ^2) );
          end
        end
      end
      
      filterKernel = filterKernel / sum(filterKernel(:));  % Normalize to sum 1
      % suppress print
      %size(filterKernel)
      
      %  PROCESS
      mat = convn( mat_in, filterKernel, 'same');
    end % gaussian3d
  end % methods (Static)
end % classdef