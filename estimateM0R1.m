
function M0R1Img = estimateM0R1( I1, I2, time1, time2, alpha1, alpha2, varargin )
  % M0R1Img = estimateM0R1( I1, I2, time1, time2, alpha1, alpha2 [ , ...
  %   'magMask', magMask, 'b1ScaleMap', b1ScaleMap ] )

  p = inputParser;
  p.addParameter( 'magMask', [], @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'b1ScaleMap', [], @isnumeric );
  p.parse( varargin{:} );
  magMask = p.Results.magMask;
  b1ScaleMap = p.Results.b1ScaleMap;

  if numel( b1ScaleMap ) == 0
    b1ScaleMap = ones( size(I1) );
  end

  Mz1 = abs(I1) ./ sin( b1ScaleMap .* alpha1 );
  Mz2 = abs(I2) ./ sin( b1ScaleMap .* alpha2 );
  M0R1Img = ( Mz2 - Mz1 ) ./ ( time2 - time1 );
  M0R1Img = max( M0R1Img, 0 );

  if numel( magMask ) > 0
    M0R1Img( magMask == 0 ) = 0;
  end

end

