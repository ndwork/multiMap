
function [t1Map,m0Map] = mri_mapT1_repeatedExcitations( dataCube, sigTimes, thetas, ...
  sliceThickness, varargin )
  % [t1Map,m0Map] = mri_mapT1_repeatedExcitations( dataCube, sigTimes, thetas, ...
  %   [ , 'b1ScaleMap', 'b1ScaleMap, 'MzAt0', MzAt0, 'mask', mask, 'verbose', verbose ] )
  %
  % Inputs:
  % dataCube - a 3D array of size MxNxK.  Each k index corresponds to a different time
  % TEs - a 1D array of size K specifying the times of each image
  %
  % Optional Inputs:
  % MzAt0 - the intial condition of the sequence
  % mask - a 2D array of size MxN.  Only map pixels with nonzero mask values.
  % verbose - scalar; info statements made if verbose is nonzero
  %
  % Outputs:
  % t1Map - a 2D array of size MxN; units are same as sigTimes
  % m0Map - a 2D array of size MxN; values are proportional to M0
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'b1ScaleMap', [], @isnumeric );
  p.addParameter( 'MzAt0', [], @isnumeric );
  p.addParameter( 'mask', [], @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'verbose', 0, @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  b1ScaleMap = p.Results.b1ScaleMap;
  MzAt0 = p.Results.MzAt0;
  mask = p.Results.mask;
  verbose = p.Results.verbose;
  %useSliceProfiles = false;
  useSliceProfiles = true;

  sData = size( dataCube );
  if numel( mask ) == 0, mask=ones( sData(1:2) ); end
  sigTimes = sigTimes(:);

  if numel( MzAt0 ) == 0, MzAt0 = zeros( sData(1:2) ); end
  if numel( MzAt0 ) == 1, MzAt0 = MzAt0 * ones( sData(1:2) ); end
  if numel( b1ScaleMap ) == 0, b1ScaleMap = ones( sData(1:2) ); end

  fminconOptions = optimoptions( 'fmincon', 'Display', 'off' );

  m0MapCols = cell( 1, sData(2) );
  t1MapCols = cell( 1, sData(2) );

  p = parforProgress( sData(2) );
  parfor i=1:sData(2)
    %if verbose ~= 0, p.progress( i, 10 ); end    %#ok<PFBNS>
    if verbose ~= 0, p.progress( i, 1 ); end    %#ok<PFBNS>

    m0MapCol = zeros( sData(1), 1 );   %#ok<PFBNS>
    t1MapCol = zeros( sData(1), 1 );

    if max( mask(:,i) ) == 0
      m0MapCols{i} = m0MapCol;
      t1MapCols{i} = t1MapCol;
      continue;
    end

    for j=1:sData(1)
      if mask(j,i)==0 || b1ScaleMap(j,i)==0, continue; end   %#ok<PFBNS>

      thisData = abs( squeeze( dataCube(j,i,:) ) );
      theseThetas = b1ScaleMap(j,i) * thetas;
      theseThetas = theseThetas(:);

      x0 = [ max( thisData ); 1000; ];
      xLB = [ 0.9 * max( thisData ./ sin(thetas(:)) ); 1; ];
      if useSliceProfiles == false
        [params,fval,fminconFlag] = fmincon( @(tmp) norm( thisData - ...
          signalModel( tmp(1), tmp(2), sigTimes, theseThetas, MzAt0(j,i) ) ), ...
          x0, [], [], [], [], xLB, [], [], fminconOptions );   %#ok<ASGLU>
      else
        [params,fval,fminconFlag] = fmincon( @(tmp) norm( thisData - ...
          signalModel_wSliceProfiles( tmp(1), tmp(2), sigTimes, theseThetas, ...
            MzAt0(j,i), sliceThickness ) ), ...
          x0, [], [], [], [], xLB, [], [], fminconOptions );   %#ok<ASGLU>
      end

      m0MapCol(j) = params(1);
      t1MapCol(j) = params(2);

      %figure; plotnice( [0; sigTimes], [MzAt0(j,i); thisData] );  
      %hold on;  plotnice( sigTimes, signalModel( ...
      %  params(1), params(2), sigTimes, theseThetas, MzAt0(j,i) ) );
      %hold on;  plotnice( sigTimes, signalModel_wSliceProfiles( ...
      %  params(1), params(2), sigTimes, theseThetas, MzAt0(j,i) ) );
      %legendnice( 'data', 'model', 'sliceProfiles' );

      %ts = [0; sigTimes];
      %figure; plotnice( ts, [ MzAt0(j,i); thisData(:) ./ sin(thetas(:)) ] );
      %M0 = params(1);  T1 = params(2);
      %Mzs = M0 - ( M0 - MzAt0(j,i) ) * exp( -ts / T1 );
      %hold on;  plotnice( ts, Mzs );
      %legendnice( 'data', 'model' );
    end

    m0MapCols{i} = m0MapCol;
    t1MapCols{i} = t1MapCol;
  end
  p.clean;

  m0Map = cell2mat( m0MapCols );
  t1Map = cell2mat( t1MapCols );
end


function sigs = signalModel( M0, T1, ts, thetas, Mz0 )
  % T1 value
  % TEs - 1D array of times
  % Mz0 - initial condition

  sigs = zeros( numel( ts ), 1 );
  Mz = Mz0;

  Mz = M0 - ( M0 - Mz ) * exp( -ts(1) / T1 );
  sigs( 1 ) = Mz * sin( thetas(1) );
  Mz = Mz * cos( thetas(1) );

  for sigIndx = 2 : numel( ts )
    Mz = M0 - ( M0 - Mz ) * exp( -( ts(sigIndx) - ts(sigIndx-1) ) / T1 );
    sigs( sigIndx ) = Mz * sin( thetas( sigIndx ) );
    Mz = Mz * cos( thetas( sigIndx ) );
  end

end


function sigs = signalModel_wSliceProfiles( M0, T1, ts, thetas, Mz0, sliceThickness )
  % M0 - (proportional to) the proton density
  % ts - 1D array of times when signal is desired
  % T1 - exponential recovery rate
  % thetas - 1D array of flip angles
  % Mz0 is the value of Mz at time 0

  tbw = 8;
  rfDuration = 3.2;  % ms

  nRF = rfDuration / 0.01;  % 4 us increments
  rf = real( dzrf( nRF, tbw, 'ex', 'ms' ) )';  % Hamming windowed sinc function
  rf = rf / sum( rf );

  rfMs_1 = make_rfMs( rf, thetas(1), tbw, rfDuration, sliceThickness );
  nVoxLayers = size( rfMs_1, 3 );
  rfMs = zeros( 3, 3, nVoxLayers, numel(thetas) );
  rfMs(:,:,:,1) = rfMs_1;
  for thetaIndx = 2 : numel( thetas )
    rfMs(:,:,:,thetaIndx) = make_rfMs( rf, thetas(thetaIndx), tbw, rfDuration, sliceThickness );
  end

  nVoxLayers = size( rfMs, 3 );
  sigs = zeros( numel( ts ), nVoxLayers );

  for layerIndx = 1 : nVoxLayers
    M = [ 0; 0; Mz0; ];

    M(3) = M0 - ( M0 - M(3) ) .* exp( -ts(1) / T1 );
    M = rfMs(:,:,layerIndx,1) * M;  % Since Mxy=0, C2R is not necessary here
    sigs( 1, layerIndx ) = M(1);  % M(1) = Mxy

    for sigIndx = 2 : numel( ts )
      M(1) = 0;  % Note: making Mxy 0 here is assuming complete T2* decay before next signal
      M(2) = 0;
      M(3) = M0 - ( M0 - M(3) ) * exp( -( ts(sigIndx) - ts(sigIndx-1) ) / T1 );

      M = rfMs(:,:,layerIndx,sigIndx) * M;
      sigs( sigIndx, layerIndx ) = M(1);  % M(1) = Mxy
    end
  end

  sigs = abs( mean( sigs, 2 ) );
end


function rfMs = make_rfMs( rf, flipAngle, tbw, rfDuration, sliceThickness )
  % Function makes RF matrices (one matrix for each slice location)
  %
  % flipAngle - flip angle in radians

  [gammaBar,gamma] = getGammaH();  % kHz/Gauss, kRad/Gauss

  nVoxLayers = round( tbw * 25 + 1 );  % Number of layers to break the voxel up into
  if mod( nVoxLayers, 2 ) == 0, nVoxLayers = nVoxLayers + 1; end

  zL = -sliceThickness * (4/2);  % Rule of thumb: simulate for 4 x slice thickness
  zH =  sliceThickness * (4/2);
  layerIndxs = -floor(nVoxLayers/2) : floor((nVoxLayers-1)/2);
  dz = ( zH - zL ) / ( nVoxLayers + mod(nVoxLayers+1,2) );
  zSlice = dz * layerIndxs;

  nRF = numel( rf );
  dtRF = rfDuration / (nRF-1);  % dtRf is time between samples of pulse
  rfExc = rf * flipAngle;
  bw = tbw / rfDuration;  % kHz
  gzrf_amplitude = bw / gammaBar / sliceThickness;  % Gauss per cm
  gz_rot = gamma * gzrf_amplitude * dtRF * ones( nRF, 1 );	% radians per cm
  splitIndx = round( numel(gz_rot) / 2 );
  gzRefocus = -gz_rot( splitIndx:end );
  gz_total = [ gz_rot; gzRefocus; ];
  rfExc_total = [ rfExc(:); zeros(size(gzRefocus)); ];

  [a,b] = abr( rfExc_total, gz_total, zSlice );
  rfMs = mri_ab2Matrix( a, b );
end

