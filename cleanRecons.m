
function cleaned = cleanRecons( recons )

  sigmaR = 0.50;
  sigmaD = 11;
  S = 31;

  [~,~,nSlices,nRecons] = size( recons );
  cleaned = zeros( size( recons ) );

  for sliceIndx = 1 : nSlices

    newReconCells = cell( nRecons, 1 );
    parfor reconIndx = 1 : nRecons

      thisRecon = squeeze( recons( :, :, sliceIndx, reconIndx ) );

      rThisRecon = real( thisRecon );
      iThisRecon = imag( thisRecon );

      newRThisRecon = bilateralFilter( rThisRecon, 'S', S, ...
       'sigmaD', sigmaD, 'sigmaR', sigmaR );
      newIThisRecon = bilateralFilter( iThisRecon, 'S', S, ...
       'sigmaD', sigmaD, 'sigmaR', sigmaR );
      
      newReconCells{ reconIndx } = newRThisRecon + 1i * newIThisRecon;
    end

    for reconIndx = 1 : nRecons
      cleaned( :, :, sliceIndx, reconIndx ) = newReconCells{ reconIndx };
    end
  end
end
