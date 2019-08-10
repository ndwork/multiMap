function run_spinEchoAnalysis
  close all; clear; rng(1);

  datacase = 8;
  xShift = -8;
  yShift = 0;
  nRot90 = 0;
  seDir = '../data/bottles4/';
  spinEchoFiles = { 'P40448.7', 'P40960.7', 'P41472.7', 'P41984.7' };
  spinEchoTimes = [ 8, 15, 30, 60 ];
  showScale = 5;

  nSpinEchoFiles = numel( spinEchoFiles );
  spinEchoData = zeros( 128, 128, nSpinEchoFiles );
  for i=1:nSpinEchoFiles
    spinEchoFile = [ seDir, spinEchoFiles{i} ];
    [rawData,header] = rawloadX( spinEchoFile );                                      %#ok<ASGLU>
    thisRecon = mri_fftRecon( rawData );
    cShiftAmount = size( thisRecon, 1 ) / 2;
    thisRecon = circshift( rot90( thisRecon, nRot90 ), [yShift cShiftAmount-xShift] );
    spinEchoData(:,:,i) = fliplr( thisRecon );
  end
  showImageCube( abs( spinEchoData ), showScale, 'border', 2, 'borderValue', 'max' );

  verbose = 0;

  outDir = [ 'output_', num2str(datacase,'%3.3i') ];
  mkdir( outDir );

  [recons,acqTimes,rfSatTime,rf30Time,seTimes,TR,tSat,tissueType,magThresh,b0Bound] = ...
    loadDatacase( datacase, 'showScale', showScale, 'verbose', verbose );    %#ok<ASGLU>

  magMask = makeMaskFromRecons( recons, magThresh );

  mmSpinEchoData = squeeze( recons( :, :, 1, 8:10 ) );
  mmSpinEchoTimes = acqTimes( 8 : end ) - tSat;
  mmRfTimes = seTimes - tSat;

  analyzeOnePoint = 0;
  if analyzeOnePoint ~= 0
    mipMagImg = max( abs( squeeze( recons( :, :, 1, : ) ) ), [], 3 );
    pt2Analyze = pickAPoint( mipMagImg, 5 );
    titlenice( 'Pick a Point to Analyze' );

    x = pt2Analyze(1);  y = pt2Analyze(2);
    spinEchoLine = squeeze( abs( spinEchoData( y, x, : ) ) );
    spinEchoLine = spinEchoLine / spinEchoLine(1);

    mmLine = squeeze( abs( mmSpinEchoData( y, x, : ) ) );
    mmLine = mmLine / mmLine(1);

    figure;  plotnice( spinEchoTimes, spinEchoLine );
    hold on; plotnice( mmSpinEchoTimes, mmLine );
    legendnice( 'spin echo', 'multimap' );

  else

    figure; showImageCube( abs( mmSpinEchoData ), 5 );
    titlenice( 'multiMap Spin Echo Data' );

    disp([ 'spin echo time: ', num2str(spinEchoTimes) ]);
    disp([ 'multiMap spin echo times: ', num2str(mmSpinEchoTimes) ]);

    mmT2Map = mri_mapT2( abs(mmSpinEchoData), mmSpinEchoTimes, 'mask', magMask );
    figure; imshowscale( mmT2Map, showScale, 'range', [0 30] );
    titlenice( 'multiMap T2' );  colorbarnice;

    t2Map = mri_mapT2( spinEchoData, spinEchoTimes, 'mask', magMask );
    figure; imshowscale( t2Map, showScale, 'range', [0 30] );
    titlenice( 'T2 Map' );  colorbarnice;

    %excFlipAngle = 30;
    %mmT2Map = multiMapT2( mmSpinEchoData, mmSpinEchoTimes, mmRfTimes, excFlipAngle, 'mask', magMask );
    %figure; imshowscale( mmT2Map, showScale, 'range', [0 100] );
    %titlenice( 'multiMap T2' );  colorbarnice;
  end

end
