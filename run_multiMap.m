
function run_multiMap( varargin )
  close all;  rng(1);
  if ~exist( 'varargin', 'var') || numel(varargin)==0, clear; varargin={}; end

  datacase = 1;  % bottles
  showScale = 5;
  verbose = true;
  outDir = [ './outputs/output_', num2str(datacase) ];


  %if ~exist( 'varargin', 'var'), varargin={}; end
  p = inputParser;
  p.addOptional( 'datacase', datacase, @ispositive );
  p.addParameter( 'outDir', outDir, @(x) true );
  p.addParameter( 'verbose', verbose, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  datacase = p.Results.datacase;
  outDir = p.Results.outDir;
  verbose = p.Results.verbose;

  if ~exist( outDir, 'dir' ), mkdir( outDir ); end

  [ recons, acqTimes, rfSatTime, rf4t1Times, seTimes, TR, tSat, sliceThickness, ...
    tissueType, magThresh, b0Bound, fieldStrength ] = loadDatacase( datacase );   %#ok<ASGLU>
  %recons = recons( :, :, 1, [ 1:3 6:13 16:20 ] );
  %acqTimes = acqTimes([ 1:3 6:end ]);

  %showAllRecons( recons, showScale, 'border', 2, 'borderValue', 'max', 'saveDir', outDir, ...
  %  'verbose', verbose );

  [db0Map, b1ScaleMap, m0Map, t1Map, t2Map, t2StarMap, ffMap] = multiMap( ...
    recons, acqTimes, rfSatTime, tSat, sliceThickness, seTimes, magThresh, fieldStrength, ...
    'b0Bound', b0Bound, 'showScale', showScale, 'verbose', verbose );

  figH = figure; imshowscale( db0Map * 1d6, showScale );  colorbarnice;
  titlenice( 'db0Map (uT)' );  saveas( gcf, [ outDir, '/db0Map.png' ] );  close( figH );

  figH = figure;  imshowscale( b1ScaleMap, showScale, 'range', [0.6 1.3] );  colorbarnice;
  titlenice( 'B1 Scale Map' );  saveas( gcf, [ outDir, '/b1ScaleMap.png' ] );  close( figH );

  figH = figure;  imshowscale( m0Map, showScale );  colorbarnice;
  titlenice( 'M0 Map' );  saveas( gcf, [ outDir, '/m0Map.png' ] );  close( figH );

  figH = figure;  imshowscale( m0Map, showScale, 'range', 'nice' );  colorbarnice;
  titlenice( 'M0 Map' );  saveas( gcf, [ outDir, '/m0Map_nice.png' ] );  close( figH );

  figH = figure;  imshowscale( t1Map, showScale, 'range', [0 300] );  colorbarnice;
  titlenice( 'T1 Map' );  saveas( gcf, [ outDir, '/t1Map.png' ] );  close( figH );

  figH = figure;  imshowscale( t1Map, showScale, 'range', [0 500] );  colorbarnice;
  titlenice( 'T1 Map' );  saveas( gcf, [ outDir, '/t1Map_largeRange.png' ] );  close( figH );

  figH = figure;  imshowscale( t1Map, showScale, 'range', 'nice' );  colorbarnice;
  titlenice( 'T1 Map' );  saveas( gcf, [ outDir, '/t1Map_nice.png' ] );  close( figH );

  figH = figure;  imshowscale( t2Map, showScale, 'range', [0 50] );  colorbarnice;
  titlenice( 'T2 Map' );  saveas( gcf, [ outDir, '/t2Map.png' ] );  close( figH );

  figH = figure;  imshowscale( t2Map, showScale, 'range', [0 120] );  colorbarnice;
  titlenice( 'T2 Map' );  saveas( gcf, [ outDir, '/t2Map_LargeRange.png' ] );  close( figH );

  figH = figure;  imshowscale( t2Map, showScale, 'range', 'nice' );  colorbarnice;
  titlenice( 'T2 Map' );  saveas( gcf, [ outDir, '/t2Map_nice.png' ] );  close( figH );

  figH = figure;  imshowscale( t2StarMap, showScale, 'range', [0 30] );  colorbarnice;
  titlenice( 'T2^* Map' );  saveas( gcf, [ outDir, '/t2StarMap.png' ] );  close( figH );

  figH = figure;  imshowscale( ffMap, showScale, 'range', [0 1] );  colorbarnice;
  titlenice( 'Fat Fraction Map' );  saveas( gcf, [ outDir, '/ffMap.png' ] );  close( figH );

  figH = figure;  imshowscale( ffMap, showScale, 'range', [0 0.4] );  colorbarnice;
  titlenice( 'Fat Fraction Map' );  saveas( gcf, [ outDir, '/ffMap_smallRange.png' ] );  close( figH );

end
