
function T1OverM0Img = estimateT1OverM0( I1, I2, ...
  time1, time2, alpha1, alpha2, magMask )


  Mz1 = abs(I1) ./ sin( alpha1 );
  Mz2 = abs(I2) ./ sin( alpha2 );
  T1OverM0Img = ( time2 - time1 ) ./ ( Mz2 - Mz1 );
  T1OverM0Img( magMask == 0 ) = 0;

end

