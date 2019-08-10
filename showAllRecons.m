
function showAllRecons( recons, varargin )

  p = inputParser;
  p.addOptional( 'showScale', 1, @isnumeric );
  p.addParameter( 'saveDir', [], @(x) true );
  p.addParameter( 'border', 0, @isnumeric );
  p.addParameter( 'borderValue', 0, @(x) true );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  showScale = p.Results.showScale;
  saveDir = p.Results.saveDir;
  border = p.Results.border;
  borderValue = p.Results.borderValue;
  verbose = p.Results.verbose;

  [Ny,Nx,nSlices,nPhases] = size( recons );  %#ok<ASGLU>
  nOutRows = 2;
  nOutCols = nPhases / nOutRows;

  figH = figure;
  showImageCube( abs( squeeze( recons(:,:,1,:) ) ), showScale, 'range', [0 1], ...
    'border', border, 'borderValue', borderValue, 'nImgsPerRow', nOutCols );
  titlenice( 'Magnitude Recons' );
  if numel( saveDir ) > 0
    saveas( gcf, [ saveDir, '/reconsMag.png' ] );
  end
  if verbose ~= true, close( figH ); end

  figH = figure;
  showImageCube( 20*log10( abs( squeeze( recons(:,:,1,:) ) ) + 1d-1 ), showScale, ...
    'border', border, 'borderValue', borderValue, 'range', 'nice', ...
    'nImgsPerRow', nOutCols );
  titlenice( 'Magnitude Recons (dB)' );
  if numel( saveDir ) > 0
    saveas( gcf, [ saveDir, '/reconsMag_dB.png' ] );
  end
  if verbose ~= true, close( figH ); end

  figH = figure;
  showImageCube( angle( squeeze( recons(:,:,1,:) ) ), showScale, ...
    'border', border, 'borderValue', borderValue, 'nImgsPerRow', nOutCols );
  titlenice( 'Phase Recons' );
  if numel( saveDir ) > 0
    saveas( gcf, [ saveDir, '/reconsPhase.png' ] );
  end
  if verbose ~= true, close( figH ); end

  figH = figure;
  showImageCube( real( squeeze( recons(:,:,1,:) ) ), showScale, ...
    'border', border, 'borderValue', borderValue, 'nImgsPerRow', nOutCols );
  titlenice( 'Real Recons' );
  if numel( saveDir ) > 0
    saveas( gcf, [ saveDir, '/reconsReal.png' ] );
  end
  if verbose ~= true, close( figH ); end

  figH = figure;
  showImageCube( imag( squeeze( recons(:,:,1,:) ) ), showScale, ...
    'border', border, 'borderValue', borderValue, 'nImgsPerRow', nOutCols );
  titlenice( 'Imag Recons' );
  if numel( saveDir ) > 0
    saveas( gcf, [ saveDir, '/reconsImag.png' ] );
  end
  if verbose ~= true, close( figH ); end

end

