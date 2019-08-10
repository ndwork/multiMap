
function [ recons, acqTimes, rfSatTime, rf4t1Times, seTimes, TR, tSat, sliceThickness, ...
  tissueType, magThresh, b0Bound, fieldStrength ] = loadDatacase( datacase, varargin )

  p = inputParser;
  p.addParameter( 'showScale', 1, @isnumeric );
  p.addParameter( 'verbose', 0, @isnumeric );
  p.parse( varargin{:} );
  showScale = p.Results.showScale;
  verbose = p.Results.verbose;

  xShift = 0;
  yShift = 0;
  b0Bound = 0.5;
  fieldStrength = 1.5;  % Tesla

  dataDir = '../data/';

  switch datacase
    case -1
      % bottles, 128 x 128, FOV 19 x 19, slice = 4.0 mm
      % optr = 2000, tSat = 1000
      tissueType = 'bottles7';
      pFile = [ dataDir, 'P27136.7' ];
      acqTimes = [ 5340 8092 10844 13596 16348 103352 303352 ...
                   1005340 1015820 1029740 1043660 ] / 1000;  % ms
      rfSatTime = 1900 / 1000;  % ms
      rf4t1Times = [ 100000, 300000 ] / 1000;  % ms
      seTimes = [ 1008860 1022780 1036700 ] / 1000;  % ms
      TR = 2000;  % ms
      tSat = 1000;  % ms
      sliceThickness = 4;  % mm
      nRot90 = 0;
      magThresh = 0.25;
      xShift = -6;
      b0Bound = 0.8;
      
    case 1
      % bottles, 128 x 128, FOV 19 x 19, slice = 4.0 mm
      % optr = 2000, tSat = 1000
      tissueType = 'bottles8';
      pFile = [ dataDir, 'P31232.7' ];
      acqTimes = [ 5324 8004 10684 13364 16044 43312 103312 1205324 ...
        1215788 1229676 1243564 ] / 1000;  % ms
      rfSatTime = 1900 / 1000;  % ms
      rf4t1Times = [ 40000, 100000 ] / 1000;  % ms
      seTimes = [ 1208844 1222732 1236620 ] / 1000;  % ms
      TR = 2400;  % ms
      tSat = 1200;  % ms
      sliceThickness = 0.4;  % cm
      nRot90 = 0;
      magThresh = 0.25;
      xShift = -7;
      b0Bound = 0.8;
  end

  nSeg = 2;
  [rawData,header] = rawloadX( pFile );   %#ok<ASGLU>
  nAcqsPerSeg = size( rawData, 4 ) / nSeg;

  if ~exist( 'shiftIndxs', 'var' ), shiftIndxs = 1:size(rawData,4); end

  % Flip the acquisitions due to spin echoes
  flipIndxs = [ nAcqsPerSeg-2 nAcqsPerSeg 2*nAcqsPerSeg-2 2*nAcqsPerSeg ];
  for sliceIndx = 1 : size( rawData, 3 )
    for dataIndx = flipIndxs
      rawData(:,:,sliceIndx,dataIndx) = fliplr( rawData(:,:,sliceIndx,dataIndx) );
    end
  end

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
    cShiftAmount = size( recons, 2 ) / 2;
    recons = circshift( recons, [0 cShiftAmount] );
  end

  recons = conj( recons );

  if xShift ~= 0 || yShift ~= 0
    for slice = 1 : size( recons, 3 )
      for reconIndx = shiftIndxs
        tmp = circshift( squeeze(recons(:,:,slice,reconIndx)), [yShift xShift] );
        recons(:,:,slice,reconIndx) = tmp;
      end
    end
  end

%   % Flip the acquisitions due to spin echoes
%   flipIndxs = [ nAcqsPerSeg-2 nAcqsPerSeg 2*nAcqsPerSeg-2 2*nAcqsPerSeg ];
%   for dataIndx = flipIndxs
%    recons(:,:,1,dataIndx) = fliplr( recons(:,:,1,dataIndx) );
%   end

  if nRot90 ~= 0
    recons = rot90( recons, nRot90 );
  end

  if verbose ~= 0
    nOutRows = 2;
    nOutCols = nPhases / nOutRows;
    showImageCube( abs( squeeze(rawData(:,:,1,:)) ), showScale, ...
      'border', 2, 'borderValue', 'max', 'nImgsPerRow', nOutCols );
    titlenice( 'Magnitude Raw Data' );
    showAllRecons( recons, showScale, 'border', 2, 'borderValue', 'max' );
  end

end

