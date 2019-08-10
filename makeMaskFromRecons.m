
function mask = makeMaskFromRecons( recons, magThresh )

  avgData = sum( recons, 3 ) / size( recons, 3 );
  mask = abs( avgData ) > magThresh;

  mask = imclose( mask, ones(3) );
  %mask = imopen( mask, ones(3) );
  %mask = imclose( mask, ones(3) );

  strelShape = strel( 'square', 3 );
  mask = imerode( mask, strelShape );

  edgeMask = makeEdgeMask( size( mask ), 2 );
  mask( edgeMask == 1 ) = 0;
end
