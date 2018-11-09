
function [ recons, acqTimes, rfSatTime, rf30Time, seTimes, TR, tSat ] = ...
  loadDatacase( datacase, varargin )

  p = inputParser;
  p.addParameter( 'showScale', 1, @isnumeric );
  p.addParameter( 'verbose', 1, @isnumeric );
  p.parse( varargin{:} );
  showScale = p.Results.showScale;
  verbose = p.Results.verbose;

  switch datacase

    case 1
      % bottles
      pFile = '../data/P50688.7' ;
      acqTimes = [ 5148 7596 10044 12492 19308 ...
                   205148 215996 228492] / 1000;  % ms
      rfSatTime = 1900 / 1000;  % ms
      rf30Time = 16060 / 1000;  % ms
      seTimes = [ 209748, 222244 ] / 1000; % ms
      TR = 400;  %ms
      tSat = 200;  % Saturation Recovery Time in ms
      nRot90 = 0;

    case 2
      % knee
      pFile = '../data/P09728.7' ;
      acqTimes = [ 5148 7596 10044 12492 19308 ...
                   205148 215996 228492] / 1000;  % ms
      rfSatTime = 1900 / 1000;  % ms
      rf30Time = 16060 / 1000;  % ms
      seTimes = [ 209748, 222244 ] / 1000; % ms
      TR = 400;  %ms
      tSat = 200;  % Saturation Recovery Time in ms
      nRot90 = 0;

    case 3
      % brain
      pFile = '../data/P17920.7' ;
      acqTimes = [ 5148 7596 10044 12492 19308 ...
                   205148 215996 228492] / 1000;  % ms
      rfSatTime = 1900 / 1000;  % ms
      rf30Time = 16060 / 1000;  % ms
      seTimes = [ 209748, 222244 ] / 1000; % ms
      TR = 400;  %ms
      tSat = 200;  % Saturation Recovery Time in ms
      nRot90 = 0;
      
    otherwise
      error('Unrecognized data case');
  end

  nSeg = 2;
  [rawData,header] = rawloadX( pFile );                                                    %#ok<ASGLU>
  [Ny,Nx,nSlices,nPhases] = size( rawData );
  nAcqsPerSeg = nPhases / nSeg;

%   % Preprocess raw data to eliminate phase offset during first aquisitions
%   for seg=1:2
%     for acq=2:4
%       thisData = rawData(:,:,1,acq+(seg-1)*nAcqsPerSeg);
%       absThisData = abs( thisData );
%       [~,maxIndx] = max( absThisData(:) );
%       [maxY,maxX] = ind2sub( size(thisData), maxIndx );
%         % Nick, whether X or Y shift depends on freq encode direction
% 
%       sData = size( absThisData );
%       %xShift = sData(2)/2+1-maxX;
%       yShift = sData(1)/2+1-maxY;
%       thisData = circshift( thisData, [yShift 0] );
%       %thisData = circshift( thisData, [0 xShift] );
%       rawData(:,:,1,acq+(seg-1)*nAcqsPerSeg) = thisData;
%     end
%   end

  % Flip the 7th acquisition due to trajectory direction
  rawData(:,:,1,7) = fliplr( rawData(:,:,1,7) );
  rawData(:,:,1,nAcqsPerSeg+7) = fliplr( rawData(:,:,1,nAcqsPerSeg+7) );

  nDataDims = ndims( rawData );
  if nDataDims == 5
    % Multiple coils
    maps = mri_makeSenseMaps( squeeze(rawData(:,:,:,1,:)) );
    sRawData = size( rawData );
    recons = zeros( sRawData(1:4) );
    for ph=1:size(rawData,4)
      phRecons = mri_fftRecon( squeeze( rawData(:,:,:,ph,:) ) );
      recons(:,:,:,ph) = rot90( mri_senseRecon( phRecons, maps ), -1 );
    end

  else
    % Single coil
    recons = mri_fftRecon( rawData );
    cShiftAmount = size( recons, 1 ) / 2;
    recons = circshift( rot90( recons, nRot90 ), [0 cShiftAmount] );
  end

  if verbose ~= 0
    nOutRows = 2;
    nOutCols = nPhases / nOutRows;
    showImageCube( abs( squeeze(rawData(:,:,1,:)) ), showScale, 'nImgsPerRow', nOutCols );
    titlenice( 'Magnitude Raw Data' );
    showImageCube( abs( squeeze( recons(:,:,1,:) ) ), showScale, 'nImgsPerRow', nOutCols );
    titlenice( 'Magnitude Recons' );
    showImageCube( angle( squeeze( recons(:,:,1,:) ) ), showScale, 'nImgsPerRow', nOutCols );
    titlenice( 'Phase Recons' );
    showImageCube( real( squeeze( recons(:,:,1,:) ) ), showScale, 'nImgsPerRow', nOutCols );
    titlenice( 'Real Recons' );
    showImageCube( imag( squeeze( recons(:,:,1,:) ) ), showScale, 'nImgsPerRow', nOutCols );
    titlenice( 'Imag Recons' );
  end
end

