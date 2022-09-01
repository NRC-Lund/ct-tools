classdef Wire < handle
  %for the moment meant to represent a single electrode wire
  properties
    manual_init;
    pp;
    I;
    tip_coords;
    wire_label;
    ROI_index;
    struct_name;
    stereo_coords;
  end %properties


  methods
    function obj = Wire( wire_label, ROI_index, coordinates, image_data)
      disp( ['Initializing wire ' ...
               num2str(wire_label) ' from coordinates']);
      % separate to check
      x = coordinates(:, 1);
      y = coordinates(:, 2);
      z = coordinates(:, 3);

      % make sure that z is decreasing in steps of 1
      if z(end) > z(1)
        x = fliplr(x);
        y = fliplr(y);
        z = fliplr(z);
      end
      assert( isequal( z', max(z)-(0:(numel(z)-1)) ) )
      
      % set properties
      obj.manual_init =  [x y z];
      obj.pp          =  [x y z];
      obj.I           = arrayfun( @( a,b,c) image_data( a,b,c), x, y, z);
      obj.wire_label  = wire_label;
      obj.ROI_index   = ROI_index;
      obj.struct_name = 'empty';
    end %Wire


    function failure = extend( obj, data)
      import electroslice.*;
      
      disp([newline 'Segmenting wire ' num2str(obj.wire_label)]);
      % watch for failure
      failure = true;

      % search parameters
      limit = size( data, 3); % z dimension, putative limit
      radius = 1; % probably doesn't make sense to change this one
      % we need more stiffness in the curves
      min_previousN = 6;
      max_previousN = 60;
      step_N = -2;
      % not too strict or it ends prematurely
      %min_condition = 0.75;
      %step = 0.02;
      %max_condition = 0.85;
      condition = 0.8; % not sure I want to change this for now
      last_pick = 1; % just to break ties, almost never used

      loop_I = obj.I( 1: size(obj.manual_init, 1));
      loop_pp = obj.manual_init;
      while size( loop_I, 1) < limit
        % local search loop
        possible_np = [];
        possible_I  = [];
        
        % loop radius of possible curves, start from stiffest!
        for tgt_previousN = max_previousN : step_N : min_previousN 
          % allowing for more previous points than manually assigned
          previousN = min( tgt_previousN, size( loop_I, 1) );
          
          % get possible np
          np_float = Wire.next_point( loop_pp, previousN);
          np = round(np_float); % round here, not fix
          
          % find maximum in vicinity
          % radius is integer and usually in 1:2
          npx = max(1, np(1)-radius) : max(1, np(1)+radius); % fix for unusual situations
          npy = max(1, np(2)-radius) : max(1, np(2)+radius); % same
          npz = np(3); % same slice

          if npz < 1
            % no possible continuation; wire has ended;
            % segmentation likely wrong but better to fix manually
            % than not to complete
            failure = false;
            break
          end
          
          slice = data( npx, npy, npz);
          M = max( slice(:));
          pos = find( slice == M);
          
          if numel( pos) == 1
            % only one option
            [a, b] = ind2sub( size( slice), pos);
            npx = npx( a);
            npy = npy( b);  
          else
            % find points among intensity M with minimal distance to 
            % (unrounded) center, will typically be the center
            [a, b] = ind2sub( size( slice), pos);
            deltas = [ a-np_float(1); b-np_float(2)];
            deltas = vecnorm( deltas);
            mm = min( deltas);
            % might still be more than one, we define a way
            % to pick while keeping deterministic behavior;
            % this will alternately pick towards bottom/left, right/top
            % ties are extremely rare anyway
            mm_ind = find( deltas == mm);
            % alternate end and 1
            last_pick = last_pick * -1;
            mm_ind = mm_ind( end*(1+last_pick)/2 + (1-last_pick)/2);

            npx = npx( a(mm_ind));
            npy = npy( b(mm_ind));
          end % numel == 1

          % compare intensities
          frac = M / Wire.runningI( loop_I, previousN );
          if frac < condition
            continue % do not keep this np
          end % if frac
          
          % check if this point is new
          distances = pdist( [possible_np; [npx npy npz]] );
          if min(distances) == 0 % also fails when distances is empty
            % it is not a new point
            continue % do not keep this np
          end % if min
          
          possible_np = [possible_np; [npx npy npz]];
          possible_I  = [possible_I; M];
        end % for tgt_previousN

        % figure out how many points were found  
        if isempty( possible_I)
          % no possible continuation; wire has ended
          failure = false;
          break
        else
          % there is at least one possible continuation
          % use the first one which is also the stiffest
          loop_I  = [loop_I;  possible_I(1)   ];
          loop_pp = [loop_pp; possible_np(1,:)];
          
          % checking for discontinuities in z
          if loop_pp(end, 3) ~= loop_pp(end-1, 3) -1
            keyboard
          end
        end % if isempty
      end % while size 

      if failure
        disp(['Attempt to segment wire ' num2str(obj.wire_label) ...
                ' failed! Proceed manually.']);
        return
      else
        disp(['Wire ' num2str(obj.wire_label) ...
                ' segmented successfully through ' ...
                num2str( size( loop_pp, 1)) ' frames!' ]);
      end % if failure
                
      % actual assignments
      obj.I  = loop_I;
      obj.pp = loop_pp;
      obj.tip_coords = Wire.get_tip_cube( loop_pp(end,:));     
      
      % ready to return
    end %extend
    
    
    function tip = get_tip( obj)
      tip = obj.tip_coords(14, :);  % center of cube
    end % get_tip
  end %methods
  
  
  methods (Static) 
    function cube = get_tip_cube( coords)
      % there is a more elegant way to write this
      % but I don't want to find it right now
      npx = coords( 1);
      npy = coords( 2);
      npz = coords( 3);

      cube = [];
      for a = (npx-1):(npx+1)
        for b = (npy-1):(npy+1)
          for c = (npz-1):(npz+1)
            cube = [cube; [a b c]];
          end
        end
      end % for a
    end % get_tip_cube
    
    
    function [centroid, out_line] = best_line( points, lastN, opt_debug)
      % get last lastN points from the current wire
      points = points( (end-lastN+1):end, : )'; % has to be row
      
      % find centroid
      centroid = mean( points, 2);
      A = points - centroid;
      [U, ~, ~] = svd(A); 
      % the left singular vectors of A are a set of
      % orthonormal eigenvectors of A*A'
      % the first one is the best fit line
      out_line = U(:,1);

      % debug plot
      if exist('opt_debug', 'var') && opt_debug == true      
        % t parametrizes the line
        % but can miss the exact slices due to rounding error
        t = d' * A;

        plot_line = centroid + [max(t) min(t)] .* d;

        hold on;
        figure;
        plot3( points(1, :), points(2, :), points(3, :), 'o');
        xlim( [min(points(1,:))-2 max(points(1,:))+2]);
        ylim( [min(points(2,:))-2 max(points(2,:))+2]);
        zlim( [min(points(3,:))-2 max(points(3,:))+2]);
        line( [plot_line(1,2), plot_line(1,1)],...
              [plot_line(2,2), plot_line(2,1)],...
              [plot_line(3,2), plot_line(3,1)], 'Color', 'black');
        grid on; 
        hold off;
      end % debug
    end %best_line


    function np = next_point( points, lastN)
      import electroslice.*;
      
      assert( size( points, 1) >= lastN)
      
      [centroid, line] = Wire.best_line( points, lastN);
      
      % get parameter that lands in next slice
      % which is always descending in z
      t = ( (points(end, 3)-1) - centroid(3)) / line(3);
      
      np = centroid + t * line;
      % notice no rounding happens here
    end %next_point
    
    
    function intensity = runningI( intensities, lastN)
      assert( numel(intensities) >= lastN)
      intensity =  mean( intensities( end-lastN+1 : end) );
    end % runningI
    
  end % methods (Static)
end %classdef
