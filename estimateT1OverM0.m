
function T1OverM0Img = estimateT1OverM0( I1, I2, time1, time2, alpha1, alpha2, varargin )
  % T1OverM0Img = estimateT1OverM0( I1, I2, time1, time2, alpha1, alpha2 [, ...
  %   'b1ScaleMap', b1ScaleMap, 'magMask', magMask ] )

  p = inputParser;
  p.addParameter( 'magMask', [], @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'b1ScaleMap', [], @isnumeric );
  p.parse( varargin{:} );
  magMask = p.Results.magMask;
  b1ScaleMap = p.Results.b1ScaleMap;

  if numel( b1ScaleMap ) == 0, b1ScaleMap = ones( size(I1) ); end

  Mz1 = abs(I1) ./ sin( b1ScaleMap .* alpha1 );
  Mz2 = abs(I2) ./ sin( b1ScaleMap .* alpha2 );
  T1OverM0Img = ( time2 - time1 ) ./ ( Mz2 - Mz1 );

  T1OverM0Img = max( T1OverM0Img, 0 );

  if numel( magMask ) > 0, T1OverM0Img( magMask == 0 ) = 0; end
end
