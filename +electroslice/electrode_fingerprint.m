%% Luciano Censoni, 2019
%% supporting routine: read spreadsheet and attempt to plot electrode fingerprint
function success = electrode_fingerprint( filename)
  %% read file
  success = false;
  try
    [num, txt, ~] = xlsread( filename, 1);
    label = [];
    ml = []; % x
    ap = []; % y
    color = [];

    start = 0;
    for k=1:size( txt, 1)
      if strcmp( txt(k, 1), '<MATLAB>TABLE STARTS HERE')
        start = k + 1;
        break
      end
    end

    if start == 0
      disp('Problem: spreadsheet has unexpected format!');
      return
    end
  catch
    disp('Problem loading or reading the file!');
    return
  end % try

  %% translate
  for k=1:size( num, 1)
    % changes to deal with AI
    if size( txt, 2) == 10
      % new xls format
      if isempty( str2num( txt{k+start, 7}))
        assert( isempty( str2num( txt{k+start, 8})) )
        continue
      end
      ap(end+1) = str2num( txt{k+start, 7});
      ml(end+1) = str2num( txt{k+start, 8});
    else
      ap(end+1) = num(k, 5);
      ml(end+1) = num(k, 6);
    end
    
    % choose color for plotting
    ventral = strcmp( char( txt(start+k, 2)), 'Ventral');
    if ventral % ventral is lateral is red
      color(end+1) = 'r';
    else % dorsal is medial is blue
      color(end+1) = 'b';
    end
    label(end+1) = num(k, 1);
  end
  
  %% now plot
  figure('Renderer', 'painters', 'Position', [0 0 1100 800]);

  hold on
  title(['Electrode ' filename ' fingerprint'], 'Interpreter', 'none');

  xlim( [ min(ml)-0.5, max(ml)+0.5]);  
  ylim( [ min(ap)-0.5, max(ap)+0.5]);
  for k = 1:numel( ap)
    plot( ml(k),ap(k),  [color(k) '.'], 'MarkerSize', 15);
    text( ml(k)+0.05, ap(k)+0.05, num2str( label(k)), 'FontSize', 7);
  end 
  hold off
  
  success = true;
  return
end
