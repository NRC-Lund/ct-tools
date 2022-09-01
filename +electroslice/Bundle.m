classdef Bundle < handle
  properties
    bundle_index;
    wireN;
    wires;
    align_helper;
  end %properties
  methods
    function obj = Bundle( wireN, planesN, bundle_index, ...
                             ROI_data, image_data, truncate, truncate_lim, min_dist)
      import electroslice.*;
      % optional adjustments
      if exist('truncate', 'var') && truncate
        image_dyn_range = double(max(image_data(:))) - double(min(image_data(:)));
        lower_lim = min(image_data(:)) + (truncate_lim(1) * image_dyn_range);
        upper_lim = min(image_data(:)) + (truncate_lim(2) * image_dyn_range);
        disp('Adjusting image intensity...');
        image_data( image_data <= lower_lim) = 0; % or lower_lim;
        image_data( image_data >= upper_lim) = upper_lim;
        disp('Done.');
      end

      if ~exist('min_dist', 'var')
        % parameters known to more or less work
        min_dist = 5; % euclidean, in pixels
      end

      max_label_dist = 0.4 * min_dist; % also euclidean, in pixels
                         
      % first: find z planes in bundle
      disp( [newline 'Locating wires in ROI index ' num2str(bundle_index)])
      coords = find( ROI_data == bundle_index);
      [x, y, z] = ind2sub( size( ROI_data), coords);
      
      disp( ['Ensuring bundle contains ' num2str(planesN) ' planes'])
      planes = sort( unique( z), 'descend')';
      assert( numel(planes) == planesN) % consistency check number of planes
      
      % looping over planes
      internal_wires = {};
      for this_plane = planes
        % for each plane
        disp( [newline 'Looping plane z = ' num2str(this_plane)])
        this_order = (z == this_plane); % indexes of coords in this plane

        this_x = x( this_order);
        this_y = y( this_order);
        this_z = z( this_order);
        
        % search area
        this_coords = sub2ind( size( ROI_data), this_x, this_y, this_z);
        search_area = numel( this_coords);
        disp( ['Search area is ' num2str( search_area)]);
        
        % target intensities
        plane_max = max( image_data( this_coords));
        last_I_frac = 0.9;
               
        % sweep whole ROI
        found_wires = 0;
        Nsweeps = 0;
        dsize = size( image_data);
        while (found_wires < wireN) && (last_I_frac > 0.6)
          disp(['Attempting frac ' num2str(last_I_frac)]);
          target_I = last_I_frac * plane_max;
          this_wires = Bundle.sweep( this_coords, min_dist, ...
                                     target_I, image_data, dsize);
          found_wires = size( this_wires, 1);
          last_I_frac = last_I_frac - 0.005; % for next round
          Nsweeps = Nsweeps + 1;
        end % while
        if ~isequal(found_wires, wireN)
          disp( [newline 'Failed to sweep plane ' num2str(this_plane) ...
            ' Found ' num2str(found_wires) ' wires with frac '...
            num2str(last_I_frac)] );
            return
        end
        
        disp(['Total sweeps for plane ' num2str(this_plane) ...
                                        ': ' num2str(Nsweeps)]);
        disp(['Final accepted intensity fraction: ' num2str(last_I_frac)]);
        
        % hopefully improve the sweeps by a quick local maximization
        disp('Starting quick local maximization...');
        maxim_round = 0;
        while true
          maxim_done = true;

          for j = 1:wireN
            wr = this_wires( j, :);
            localI = image_data( wr(1), wr(2), wr(3));
            search_x = wr(1)-1 : wr(1)+1;
            search_y = wr(2)-1 : wr(2)+1;
            search_z = wr(3); % same plane

            slice = image_data( search_x, search_y, search_z);
            M = max( slice(:));
            if M <= localI
              % already at local maximum
              %disp(['Wire ' num2str( j) ' is already at local maximum.']);
              continue % nothing to do
            end

            % get actual position
            pos = find( slice == M);
            pos = pos(1); % simplest tiebreak
            [a, b] = ind2sub( [3,3], pos);
            search_x = search_x( a);
            search_y = search_y( b);

            % see if it's still valid
            temp_wires = this_wires;
            temp_wires( j, :) = [search_x, search_y, search_z];
            distances = pdist(temp_wires);
            if min(distances) > min_dist
              % change this wire
              %disp(['Moving wire ' num2str( j) ' once towards higher intensity!']);
              this_wires( j, :) = [search_x, search_y, search_z];
              maxim_done = false; % we made at least one change
              %else
              %disp(['Moving wire ' num2str( j) ' would violate minimum distance.']);
            end % if min
          end % for j

          % check if no changes were made
          if maxim_done
            disp('Local maximization finished.');
            break
          else
            maxim_round = maxim_round + 1;
            disp(['Local maximization round ' num2str( maxim_round) ' completed.']);
          end % f maxim
        end % while

        internal_wires{end+1} = this_wires;
      end % for planes
      
      % reorder wires to maintain labels
      for k = 2:planesN
        disp( [newline 'Reordering labels in plane number ' num2str(k)])

        % elegant
        [can_solve, solved_order] = Bundle.match_labels( internal_wires{k-1}, ...
                                  internal_wires{k}, max_label_dist);
        if ~can_solve
          disp(['Unable to find adequate labels for plane ' num2str( k) '. Stopping!']);
          return
          % failed to get correct labels; write error message
        else
          disp(['Correct ordering found for plane ' num2str( k) '.' ]);
          internal_wires{k} = internal_wires{k}( solved_order, :);
        end
      end % for k
      
      % updating object property
      obj.bundle_index = bundle_index;
      obj.wireN = wireN;
      disp([newline 'Building Wire objects.']); % for good measure
      obj.wires = [];
      for k = 1:wireN
        this_wire = [];
        for j = 1:planesN
          this_wire = [this_wire; internal_wires{j}(k, :)];
        end
        % wires is a horizontal array for ease of looping
        obj.wires = [obj.wires Wire( k, bundle_index + k, ...
                                        this_wire, image_data)];
      end
      %obj.plot_labels(); % dont plot before correcting labels
      disp("Done!");
    end % Bundle
    
    
    function changed = correct_labels( obj)
      disp("Correcting wire labels.");
      % I'm imposing numerical labels for now
      changed = [];
      for wr = obj.wires
        new_label = input(['Input the new numerical label for wire ' ...
                        num2str(wr.wire_label) ', Return to skip or -1 to end: ']);
        if isempty( new_label)
          disp('Skipping!');
          continue
        elseif new_label < 0
          disp('Finished!');
          break
        end % if isempty
        changed( end+1) = wr.wire_label;
        wr.wire_label = new_label;
      end % for wr
    end % correct_labels

    
    function extend_all( obj, data, truncate, truncate_lim)
      % optional adjustments
      if exist('truncate', 'var') && truncate
        %imlook4d_current_handles = evalin('base', 'imlook4d_current_handles');
        image_dyn_range = double(max(data(:))) - double(min(data(:)));
        lower_lim = min(data(:)) + (truncate_lim(1) * image_dyn_range);
        upper_lim = min(data(:)) + (truncate_lim(2) * image_dyn_range);
        disp('Adjusting image intensity...');
        data( data <= lower_lim) = 0; % or lower_lim;
        data( data >= upper_lim) = upper_lim;
        disp('Done.');
      end

      failures = arrayfun( @(wr) wr.extend( data), obj.wires);
      
      % check success
      if sum(failures) > 0
        disp( ['Unable to extend wire(s): ' num2str( find(failures))]);
        disp( 'Set ROIs manually!');
      end
      obj.check_overlap();
    end % extend_all
    
    
    function wires = check_overlap( obj)
      % assuming wires have been extended,
      % checking quality of results
      
      % larger than this means no overlap on ROI cubes
      min_dist = sqrt(12);

      % matrix with centers of tip positions
      coords = cell2mat(arrayfun( @(n) obj.wires(n).get_tip()', ...
                          1:obj.wireN, 'UniformOutput', false))';

      % distance matrix
      distances = pdist( coords, 'euclidean');
      if min( distances) > min_dist
        wires = 0;
      else
        distances = squareform( distances);
        % just to make sure I don't count loops
        for k = 1:size( distances, 1)
          distances(k,k) = max(distances(:));
        end % for k
        
        % plot results
        figure;
        G = graph( distances <= min_dist);
        g = plot( G);
        g.NodeLabel = [obj.wires.wire_label];
        
        % return values        
        [a, b] = ind2sub( size(distances), find( distances < min_dist));
        wires = unique([a; b])'; % replace with actual wire labels
      end % if
    end % check_overlap
    
    
    function plot_labels( obj)
      %plot to verify
      figure('Renderer', 'painters', 'Position', [300 100 1000 900]);
      title(['Plane ' num2str( obj.wires(1).manual_init(end,3) ) ]);

      hold on;
      for wr = obj.wires
        this_wire = wr.manual_init(end, :); % for clarity
        plot( this_wire(1), this_wire(2), 'rO', 'MarkerSize', 3, 'MarkerFaceColor', 'r');
        text( this_wire(1)*1.002, this_wire(2)*1.002, num2str( wr.wire_label), 'FontSize', 6);
      end

      labels = 'xyz';
      xlabel( labels(1)); ylabel( labels(2));
      hold off;
    end % plot_labels


    function plot_all( obj)
      % visual aid to check for wire crossings
      figure;
      hold on;
      for wr = obj.wires
        body = wr.pp; % for clarity
        plot3( body(:, 1), body(:, 2), body(:, 3) );
      end % for wr

      labels = 'xyz';
      xlabel( labels(1)); ylabel( labels(2)); zlabel( labels(3));
      grid();
      hold off;
    end % plot_all


    function plot_tips( obj)
      % visual aid to check for wire crossings
      figure;
      hold on;
      for wr = obj.wires
        tip = wr.get_tip();
        scatter3( tip(1), tip(2), tip(3) );
      end % for wr

      labels = 'xyz';
      xlabel( labels(1)); ylabel( labels(2)); zlabel( labels(3));
      grid();
      hold off;
    end % plot_all


    function [ROI_data, ROI_names] = update_ROI( obj, ROI_data, ROI_names)
      % last_name is usually "Add ROI"
      last_name = ROI_names{ end};
      ROI_names( end) = [];
      
      % clear all ROIs from obj.bundle_index onwards
      disp('Clearing ROIs! This takes a while...');
      ROI_names = { ROI_names{ 1:obj.bundle_index } }; 
      ROI_data( ROI_data > obj.bundle_index) = 0;

      % loop wires
      for wr = obj.wires
        % update the names
        ROI_names{end+1}  = ['Electrode Tip ' num2str(wr.wire_label)];
        
        for s = 1:size( wr.tip_coords, 1)
          ROI_data( wr.tip_coords(s, 1), wr.tip_coords(s, 2), ...
                      wr.tip_coords(s, 3)) = wr.ROI_index;
        end % for s
      end % for wr
      
      % reintroduce last name
      ROI_names{end+1} = last_name;
      ROI_names = ROI_names';
      disp('Don''t forget to Import now!');
    end % update_ROI
    
    
    function changed = read_manual_ROI( obj, ROI_data, scope)
      % read manually set/corrected tip positions from ROIs
      % remember that overlaps in the ROI cubes can result in ROI CGs
      % moving slightly! A cycle of update_ROI -> read_manual_ROI ->
      % update_ROI -> read_manual_ROI -> ... is NOT equivalent to
      % identity (no operation) even if no changes are made!
      % I'm still thinking how to improve this. Meanwhile, only call
      % read_manual_ROI if actual changes have been made, and always
      % clear a given ROI before overwriting.
      
      % operation type
      import eletroslice.*;

      if strcmp( scope , 'local')
        disp('Inspecting ROIs for Local changes!');
        gl_change = false;
      elseif strcmp( scope, 'global')
        disp('Inspecting ROIs for Global changes!');
        gl_change = true;
      else
        disp('Unknown operation! Aborting...');
        changed = -1;
        return
      end % if strcmp
      
      % reading ROIs
      disp(['Reading ROIs! This takes a while...' newline]);
      all_roi = getRoiCGs( ROI_data);
      
      % results
      changed = [];
      for k = 1:obj.wireN % I will need the index later
        wr = obj.wires( k);
        
        % read tip from ROIs
        disp(['Updating wire ' num2str( wr.wire_label) ]);
        
        % rounding is necessary because of eventual ROI overlaps
        read_tip = round( all_roi( wr.ROI_index, :)); % round not fix
        
        % tip coordinates should never go from finite to NaN;
        % if it's NaN, then it was already NaN before
        % and should not change
        % exception: very dramatic rotations -- how to deal with that?
        if any( isnan( read_tip))
          disp(['Wire ' num2str( wr.wire_label) ...
                  ' is NaN and will not update now.']);
          continue
        end % if any isnan
        
        % compare
        this_tip = wr.get_tip();
        if ~all( read_tip == this_tip )
          disp([newline 'Wire ' num2str( wr.wire_label) ' was changed. Overwriting!'] );
          to_update = [ k];
          
          % if it's a global change, check who to update together before changing;
          % check only wires of lesser index and exact overlap
          if gl_change
            for nn = 1:(k-1) % lesser index
              if all( obj.wires( nn).get_tip() == this_tip) % exact overlap
                to_update( end+1) = nn;
                disp(['Updating wire ' num2str( obj.wires( nn).wire_label) ...
                  ' together with wire ' num2str( wr.wire_label) '!']);
              end % if ==
            end % for nn
          end % if gl_change
          
          for nn = to_update
            % helper function creates cube around CG instead.
            obj.wires( nn).tip_coords = Wire.get_tip_cube( read_tip);
            changed( end+1) = obj.wires( nn).wire_label;
          end % for nn
        else
          disp(['No change in wire ' num2str( wr.wire_label) '.']);
        end % if ~=
      end % for k
      
      if gl_change
        % check that all wires have been changed
        if ~all( sort( [obj.wires.wire_label]) == sort( changed) )
          disp('Attempted global procedure but could not update all wires! Check results!');
        end % if ~=
      end % if gl_change
      disp('Don''t forget to run update_ROI now!');
    end % read_manual_ROI
    
    
    function save_bundle( obj, filename)
      save( filename, 'obj', '-v7.3');
    end % save_bundle


    function create_aligner_atlas( obj, n_ROIs, structure_names, segmented_atlas, scan_resolution)
      % scan resolution is needed in both cases so it's the first check
      import electroslice.*;

      if ~exist('scan_resolution', 'var')
        scan_resolution = -1; % sentinel value
      end
      obj.align_helper = Aligner( false, scan_resolution, n_ROIs, structure_names, segmented_atlas);
    end


    function create_aligner_observed( obj, scan_resolution)
      % scan resolution is needed in both cases so it's the first check
      import electroslice.*;

      if ~exist('scan_resolution', 'var')
        scan_resolution = -1; % sentinel value
      end
      obj.align_helper = Aligner( true, scan_resolution);
    end % create_aligner_observed


    function get_names( obj)
      names = obj.align_helper.internal_get_names( obj.wires);
      assert( numel( names) == obj.wireN);
      
      % display and fill array in the same loop
      disp(['WIRE' char(9) char(9) 'STRUCTURE']);
      for k = 1:obj.wireN
        wr = obj.wires( k);
        wr.struct_name = names{ k};
        disp([num2str( wr.wire_label) char(9) char(9) wr.struct_name]);
      end
    end


    function get_stereo_coords( obj)
      coords = obj.align_helper.internal_get_stereo_coords( obj.wires);
      
      mcoords = cell2mat( coords(:));
      min_ml = min( mcoords(:, 1));
      max_ml = max( mcoords(:, 1));
      min_ap = min( mcoords(:, 2));
      max_ap = max( mcoords(:, 2));
      min_dv = min( mcoords(:, 3));
      max_dv = max( mcoords(:, 3));

      % display and fill array in the same loop
      figure;
      hold on;
      tab = char(9);
      disp(['WIRE' tab tab 'ML' tab tab tab 'AP' tab tab tab 'DV']);
      for k = 1:obj.wireN
        wr = obj.wires( k);
        wr.stereo_coords = coords{ k};
        disp([num2str( wr.wire_label) ...
                         tab tab num2str( wr.stereo_coords(1)) tab ...
                         tab tab num2str( wr.stereo_coords(2)) tab ...
                         tab tab num2str( wr.stereo_coords(3))]);
        scatter3( wr.stereo_coords(1), wr.stereo_coords(2), wr.stereo_coords(3));
      end % for
      disp(['MIN' tab tab num2str( min_ml) tab tab tab num2str( min_ap) tab tab tab num2str( min_dv)]);
      disp(['MAX' tab tab num2str( max_ml) tab tab tab num2str( max_ap) tab tab tab num2str( max_dv)]);

      % finish plot
      xlabel('ML'); ylabel('AP'); zlabel('DV')
      quiver3(0,0,0, 1, 0, 0);
      quiver3(0,0,0, 0, 1, 0);
      quiver3(0,0,0, 0, 0, 1);
      text(1.05, 0, 0, 'ML');
      text(0, 1.05, 0, 'AP');
      text(0, 0, 1.05, 'DV');
      grid();
      hold off;
    end
  end % methods


  methods (Static)
    function wires = sweep( this_coords, min_dist, target_I, ...
                                image_data, dsize)
      wires = [];
      for c = this_coords'
        % first test intensity
        if image_data( c) > target_I
          % then test collisions 
          [x, y, z] = ind2sub( dsize, c);
          temp_wires = [wires; [x, y, z]];
          distances = pdist(temp_wires);
          % decide
          if isempty(distances) || (min(distances) > min_dist)
            wires(end+1, :) = [x, y, z];
          end % if min
          
        end % if image_data
      end % for c
    end % sweep


    function [can_solve, solved_order] = match_labels( pos1, pos2, max_dist)
      can_solve = false;
      solved_order = [];

      dmat = pdist2(pos1, pos2) <= max_dist;

      if det( double(dmat)) == 0, return, end
      % simple check for linear dependence or all-zero row/column
      
      sorder = matchpairs( -1.0 * dmat, 1e9); % minimize negative
      
      % make sure the result is sorted based on the first plane
      [~, ind] = sort( sorder(:,1) );
      sorder = sorder(ind,:);

      if trace( pdist2( pos1(sorder(:,1),:), pos2(sorder(:,2),:) ) <= max_dist) == size(pos1, 1)
        can_solve = true;
        solved_order = sorder(:,2);
      end
    end

    
    function obj = load_bundle( filename)
      l = load( filename);
      obj = l.obj;
    end % load_bundle
  end % methods (Static)
end % classdef
