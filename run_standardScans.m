
function run_standardScans
  close all; clear;

  datacase = 1;
  showScale = 5;
  verbose = false;
  saveDir = './outputs/bottleTruth';

  data = loadCalibrationData( datacase );
  if ~exist( saveDir, 'dir' ), mkdir( saveDir ); end

  % data mask
  t1Data_ir = makeDataCubeFromFiles( data.t1Files_ir, data.yxShift );
  mask = makeMaskFromRecons( t1Data_ir, data.magThresh );
  if verbose ~= 0 || exist( 'saveDir', 'var' )
    figure; imshowscale( mask, showScale );  titlenice('SSMRF Mask');
  end
  if exist( 'saveDir', 'var' )
    saveas( gcf, [ saveDir, '/mask.png' ] );
    if ~verbose, close all; end
  end


%   % off resonance
%   disp([ 'NICK, YOU NEED TO MAKE SURE THAT THE CENTER FREQUENCY OF OFF RESONANCE', ...
%     ' IS THE SAME AS THAT OF THE MULTIMAP THE CENTER FREQUENCY IS IN THE PFILE' ]);
%   offResFileIndx = 1;
%   offResData = makeDataCubeFromFiles( data.offResFiles{offResFileIndx}, data.yxShift );
%   nOffRes = size( offResData, 3 );
%   offResTEs_ms = zeros( nOffRes, 1 );
%   [~,dataHeader] = rawloadX( data.offResFiles{offResFileIndx} );
%   offResTEs_ms(1) = dataHeader.te / 1000;
%   offResTEs_ms(2) = dataHeader.te2 / 1000;
%   if verbose ~= 0
%     figure;
%     showImageCube( angle(offResData), showScale, 'border', 4, 'borderValue', 'max' );
%     titlenice('offres Data (rad)');
%   end
%   %[ offResMap_radperms, phaseOffsetMap_rad ] = mri_mapOffRes( offResData, offResTEs_ms, ...
%   %  'mask', mask );
%   [ offResMap_radperms, phaseOffsetMap_rad ] = mri_mapOffResSimple( offResData, ...
%     offResTEs_ms, 'mask', mask );
%   offResMap_kHz = offResMap_radperms / (2*pi);
%   if verbose ~= 0 || exist( 'saveDir', 'var' )
%     offResFig = figure; imshowscale( mask.*offResMap_kHz, showScale );
%     titlenice('Off Resonance (kHz)');  colorbarnice;
%     figure; imshowscale( mask.*phaseOffsetMap_rad * 180/pi, showScale );
%     titlenice( 'PhaseOffset (rad)' );  colorbarnice;
%   end
%   if exist( 'saveDir', 'var' )
%     saveas( offResFig, [ saveDir, '/offResMap_kHz.png' ] );
%     if ~verbose, close all; end
%   end


  % t1 - inversion recovery
  t1Data_ir = makeDataCubeFromFiles( data.t1Files_ir, data.yxShift );
  t1MapIR_ms = mri_mapT1InversionRecovery( t1Data_ir, data.t1TIs_ir, ...
    'mask', mask, 'verbose', 1 );
  if verbose ~= 0 || exist( 'saveDir', 'var' )
    figure;
    showImageCube( abs(t1Data_ir), showScale, 'border', 4, 'borderValue', 'max' );
    titlenice('t1 ir Data');
    t1Fig = figure; imshowscale( mask.*t1MapIR_ms, showScale, 'range', [0 500] );
    titlenice( 'T1 Map IR (ms)' );  colorbarnice;
    figure; imshowscale( mask.*t1MapIR_ms, showScale, 'range', [0 5000] );
    titlenice( 'T1 Map IR (ms)' );  colorbarnice;
  end
  if exist( 'saveDir', 'var' )
    saveas( t1Fig, [ saveDir, '/t1MapIR_ms.png' ] );
    if ~verbose, close all; end
  end


  % b1
  b1Data = makeDataCubeFromFiles( data.b1Files, data.yxShift );
  divData = abs(b1Data(:,:,2)) ./ abs(b1Data(:,:,1));
  divData = max( min( divData, 2 ), 0 );
  angleMapSimple = acos( divData * 0.5 );
  b1MapSimple = angleMapSimple / min( data.b1Angles );
  if verbose ~= 0 || exist( 'saveDir', 'var' )
    figure;
    showImageCube( abs(b1Data), showScale, 'border', 4, 'borderValue', 'max' );
    titlenice('b1 Data');
    figure; imshowscale( mask.*angleMapSimple * 180/pi, showScale );
    titlenice( 'Simple B1 (deg)' );  colorbarnice;
    figure; imshowscale( mask.*b1MapSimple, showScale );
    titlenice( 'Simple B1 Map' );  colorbarnice;
  end

  b1Map = mri_doubleAngleMapB1( b1Data, data.b1Angles, 'mask', mask );
  if verbose ~= 0 || exist( 'saveDir', 'var' )
    figure; imshowscale( mask .* b1Map, showScale );
    titlenice( 'B1 Map' );  colorbarnice;
  end
  if exist( 'saveDir', 'var' )
    saveas( gcf, [ saveDir, '/b1Map.png' ] );
    if ~verbose, close all; end
  end


  % t2
  t2Data = makeDataCubeFromFiles( data.t2Files, data.yxShift );
  nT2 = size( t2Data, 3 );
  t2TEs_ms = zeros( nT2, 1 );
  for i=1:nT2
    [~,dataHeader] = rawloadX( data.t2Files{i} );
    t2TEs_ms(i) = dataHeader.te / 1000;
  end
  t2Map_ms = mri_mapT2( t2Data, t2TEs_ms, 'mask', mask, 'verbose', verbose );
  if verbose ~= 0 || exist( 'saveDir', 'var' )
    figure;
    showImageCube( abs(t2Data), showScale, 'border', 2, 'borderValue', 'max' );
    titlenice('t2 Data');
    t2Fig = figure; imshowscale( mask.*t2Map_ms, showScale, 'range', [0 50] );
    titlenice('T2 Map (ms)'); colorbarnice;
    t2Fig2 = figure; imshowscale( mask.*t2Map_ms, showScale, 'range', [0 120] );
    titlenice('T2 Map (ms)'); colorbarnice;
    t2FigNice = figure; imshowscale( mask.*t2Map_ms, showScale, 'range', 'nice' );
    titlenice('T2 Map (ms)'); colorbarnice;
  end
  if exist( 'saveDir', 'var' )
    saveas( t2Fig, [ saveDir, '/t2Map_ms.png' ] );
    saveas( t2Fig2, [ saveDir, '/t2Map_LargeRange.png' ] );
    saveas( t2FigNice, [ saveDir, '/t2FigNice.png' ] );
    if ~verbose, close all; end
  end


  % t2 Star
  t2StarData = makeDataCubeFromFiles( data.offResFiles, data.yxShift );
  t2StarTEs_ms = zeros( size( t2StarData, 3 ), 1 );
  for fileIndx=1:numel( data.offResFiles )
    [~,dataHeader] = rawloadX( data.offResFiles{fileIndx} );
    t2StarTEs_ms(2*(fileIndx-1)+1) = dataHeader.te / 1000;
    t2StarTEs_ms(2*(fileIndx-1)+2) = dataHeader.te2 / 1000;
  end
  t2StarData = cat( 3, t2StarData(:,:,1), t2StarData(:,:,2:2:end) );
  t2StarTEs_ms = [ t2StarTEs_ms(1); t2StarTEs_ms(2:2:end); ];
  t2StarMap_ms = mri_mapT2( t2StarData, t2StarTEs_ms, 'mask', mask, ...
    'verbose', verbose );
  if verbose ~= 0 || exist( 'saveDir', 'var' )
    showImageCube( abs(t2StarData), showScale, 'border', 4, 'borderValue', 'max' );
    figure; imshowscale( mask .* t2StarMap_ms, showScale, 'range', [0 50] );
    titlenice( 'T2 Star Map (ms)' );  colorbarnice;
  end
  if exist( 'saveDir', 'var' )
    saveas( gcf, [ saveDir, '/t2StarMap_ms.png' ] );
    if ~verbose, close all; end
  end

end


function data = loadCalibrationData( datacase )

  fieldStrength = 1.5;

  switch datacase
    case 1
      dataDir = '../data/bottles7/';
      offResFiles = { 'P14336.7', 'P14848.7' };
      b1Files = { 'P13312.7', 'P12800.7' };
      b1Angles = [ 60, 120 ] * pi/180;  % radians
      t2Files = { 'P16384.7', 'P16896.7', 'P17408.7', 'P17920.7', ...
        'P18432.7', 'P18944.7', 'P19456.7', 'P19968.7' };
      t1Files_ir = { 'P21504.7', 'P23552.7', 'P22016.7', ...
        'P24064.7', 'P22528.7', 'P23040.7' };
      t1TIs_ir = [ 50 80 150 300 500 1000 ];
      magThresh = 0.2;
      xShift = -6;
      b0Bound = 0.8;


    case -1
      dataDir = '../data/bottles5/';
      offResFiles = { 'P01536.7', 'P02048.7' };
      b1Files = { 'P00000.7', 'P64512.7' };
      b1Angles = [ 60, 120 ] * pi/180;  % radians
      t2Files = { 'P03072.7', 'P03584.7', 'P04096.7', 'P05120.7', ...
        'P05632.7', 'P06144.7', 'P06656.7', 'P07168.7' };
      t1Files_ir = { 'P09216.7', 'P09728.7', 'P10240.7', ...
        'P10752.7', 'P11264.7', 'P11776.7', 'P12288.7' };
      t1TIs_ir = [ 50 100 200 400 800 1200 1600];  % ms
      magThresh = 0.2;
      xShift = -7;
      b0Bound = 0.8;

    otherwise
      error('Datacase not recognized');
  end

  data.b0Bound = b0Bound;
  data.magThresh = magThresh;
  data.fieldStrength = fieldStrength;

  if ~exist( 'xShift', 'var' ), xShift=0; end
  if ~exist( 'yShift', 'var' ), yShift=0; end
  data.yxShift = [ yShift xShift ];


  for i=1:numel(offResFiles)
    offResFiles{i} = [ dataDir, offResFiles{i} ];
  end
  data.offResFiles = offResFiles;

  for i=1:numel(b1Files)
    b1Files{i} = [ dataDir, b1Files{i} ];
  end
  data.b1Files = b1Files;
  data.b1Angles = b1Angles;

  for i=1:numel(t2Files)
    t2Files{i} = [ dataDir, t2Files{i} ];
  end
  data.t2Files = t2Files;

  for i=1:numel(t1Files_ir)
    t1Files_ir{i} = [ dataDir, t1Files_ir{i} ];
  end
  data.t1Files_ir = t1Files_ir;
  data.t1TIs_ir = t1TIs_ir;  % in ms
end


function dataCube = makeDataCubeFromFiles( files, yxShift )
  % files is a 1D cell array of filenames
  if ~iscell( files ), files={files}; end

  nFiles = numel( files );
  nEchoes = 0;
  for i=1:nFiles
    [~,dataHeader] = rawloadX( files{i} );
    nEchoes = nEchoes + dataHeader.nechoes;
  end
  dataCube = zeros( dataHeader.nframes, dataHeader.frsize, nEchoes );
  shiftAmount = dataHeader.frsize / 2;
  dataIndx = 0;
  for i=1:nFiles
    if ~exist( files{i}, 'file' ), error(['File ', files{i}, ' doesn''t exist']); end
    [kData,kDataHdr] = rawloadX( files{i} );
    nEchoes = kDataHdr.nechoes;
    for j=1:nEchoes
      dataIndx = dataIndx + 1;
      dataCube(:,:,dataIndx) = mri_fftRecon( kData(:,:,j) );
      dataCube(:,:,dataIndx) = fliplr( ...
        circshift( dataCube(:,:,dataIndx), [0,shiftAmount]-yxShift ) );
    end
  end
end


