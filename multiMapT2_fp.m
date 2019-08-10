
function t2Map = multiMapT2_fp( dataCube, dataTimes, rfTimes, excFlipAngle, varargin )
  % t2Map = multiMapT2( dataCube, dataTimes, rfTimes, excFlipAngle [, b1ScaleMap, ...
  %   'mask', mask, 'verbose', verbose ] )
  %
  % Inputs:
  % excFlipAngle - the flip angle of the RF pulse in degrees
  % dataCube - an M x N x K array (where K is the number of data acquisitions)
  % dataTimes - a K element array sepcifying the times that the data was collected
  %  (for K > 1, these should be the spin echo times)
  % rfTimes - a K-1 element array specifying the times that the inversion RF pulses were 
  %   The time is the center of the rf pulse.
  %   It is assumed that the center of the excitation RF pulse is time 0
  %
  %
  % Written by Nicholas Dwork, Copyright 2019

  p = inputParser;
  p.addOptional( 'b1ScaleMap', [], @isnumeric );
  p.addParameter( 'mask', [], @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'verbose', 0, @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  b1ScaleMap = p.Results.b1ScaleMap;
  mask = p.Results.mask;
  verbose = p.Results.verbose;

  if numel( b1ScaleMap ) == 0
    b1ScaleMap = ones( size( dataCube, 1 ), size( dataCube, 2 ) );
  end
  
  tbw = 8;
  rfDuration = 3.2;  % ms
  sliceThickness = 0.5; % cm
  B0 = 0;  % residual B0 offset

  nRF = tbw * 100 + 1;
  rf = real( dzrf( nRF, tbw, 'ex', 'ms' ) )';  % Hamming windowed sinc function
  rf = rf / sum( rf );

  dtRF = rfDuration / (nRF-1);  % dtRf is time between samples of pulse
  nVoxLayers = round( tbw * 100 + 1 );  % Number of layers to break the voxel up into
  if mod( nVoxLayers, 2 ) == 0, nVoxLayers = nVoxLayers + 1; end
  zL = -sliceThickness * (4/2);  % Rule of thumb: simulate for 4 x slice thickness
  zH =  sliceThickness * (4/2);
  layerIndxs = -floor(nVoxLayers/2) : floor((nVoxLayers-1)/2);
  dz = ( zH - zL ) / ( nVoxLayers + mod(nVoxLayers+1,2) );
  zSlice = dz * layerIndxs;

  b1Scales = linspace( 0.5, 1.5, 100 );
  nB1Scales = numel( b1Scales );

  R2C = mri_makeSigConversionMatrix_r2c();
  C2R = mri_makeSigConversionMatrix_c2r();
  [gammaBar,gamma] = getGammaH();  % kHz/Gauss, kRad/Gauss

  rfExc = rf * excFlipAngle * pi/180;
  bw = tbw / rfDuration;  % kHz
  gzrf_amplitude = bw / gammaBar / sliceThickness;  % Gauss per cm
  gz_rot = gamma * gzrf_amplitude * dtRF * ones( nRF, 1 );	% radians/cm
  splitIndx = round( numel(gz_rot)/2 );
  gzRefocus = -gz_rot( splitIndx:end );
  gz_total = [ gz_rot; gzRefocus; ];
  rfExc_total = [ rfExc(:); zeros(size(gzRefocus)); ];

  rfInv = rf * pi;    % TODO: use B1 map to get more accurate flip angles here
  gzInv = gamma * gzrf_amplitude * dtRF * ones( nRF, 1 );	% radians/cm

  dRfs2TEs = dataTimes - rfTimes;

  T2s = 1:1:1000;
  nT2s = numel( T2s );
  nTEs = size( dataCube, 3 );
  simLayerSigs = zeros( nB1Scales, nT2s, nTEs, nVoxLayers );
    % B1, T2, TE, layer

  nB1Scales = numel( b1Scales );
  for b1ScaleIndx = 1 : nB1Scales
    if verbose ~= 0 && mod( b1ScaleIndx, 10 )==0
      disp([ 'Working on b1ScaleIndx ', num2str(b1ScaleIndx), ' of ', num2str( nB1Scales ) ]);
    end

    b1Scale = b1Scales( b1ScaleIndx );

    [aExc,bExc] = abr( b1Scale * rfExc_total, gz_total, zSlice, B0 );
    rfExcMs = mri_ab2Matrix( aExc, bExc );
    [aInv,bInv] = abr( b1Scale * rfInv, gzInv, zSlice, B0 );
    rfInvMs = mri_ab2Matrix( aInv, bInv );
    rfNegInvMs = mri_ab2Matrix( aInv, -bInv );

    b1SimLayerSigs = cell( 1, 1, nVoxLayers );
    parfor layer=1:nVoxLayers
      thisSimLayerSig = zeros( nT2s, nTEs );

      for t2Indx = 1 : nT2s
        T2 = T2s( t2Indx );   %#ok<PFBNS>

        M = C2R * rfExcMs(:,:,layer) * R2C * [0; 0; 1];

        for invIndx = 1:numel( dRfs2TEs )
          dRf2TE = dRfs2TEs( invIndx );        

          M(1:2) = M(1:2) * exp( -dRf2TE / T2 );

          if mod( invIndx, 2 ) == 0
            M = C2R * rfNegInvMs(:,:,layer) * R2C * M;
          else
            M = C2R * rfInvMs(:,:,layer) * R2C * M;
          end

          thisSimLayerSig( t2Indx, invIndx ) = M(1) + 1i * M(2);
        end
      end

      b1SimLayerSigs{ 1, 1, layer } = thisSimLayerSig;      
    end

    simLayerSigs(b1ScaleIndx,:,:,:) = cell2mat( b1SimLayerSigs );
  end
  simLayerSigs = single( simLayerSigs );
save( 'simLayerSigs.mat', 'simLayerSigs', '-v7.3' );
load simLayerSigs.mat;

  t2SimSigs = abs( sum( simLayerSigs, ndims(simLayerSigs) ) );
  t2Norms = norms( t2SimSigs, 2, 2 );
  t2SimSigs = t2SimSigs ./ repmat( t2Norms, [ 1 size(t2SimSigs,2) 1 ] );

  sData = size( dataCube );
  t2Map = cell( 1, sData(2) );

  p = parforProgress( sData(2) );
  parfor i = 1 : size( dataCube, 2 )
    if verbose ~= 0, p.progress( i, 10 ); end  %#ok<PFBNS>

    thisMapCol = zeros( size( dataCube, 1 ), 1 );
    for j = 1 : size( dataCube, 1 )
      if mask(j,i) == 0, continue; end

      dataVec = squeeze( abs( dataCube(j,i,:) ) );
      normdVec = dataVec ./ norm( dataVec, 2 );

      thisB1Scale = b1ScaleMap( j, i );
      [~,b1ScaleIndx] = min( abs( thisB1Scale - b1Scales ) );
      t2SimSigs4ThisB1Scale = squeeze( t2SimSigs( b1ScaleIndx, :, : ) );   %#ok<PFBNS>
      normSims = norms( t2SimSigs4ThisB1Scale, 2, 2 );
      normdSimSigs = t2SimSigs4ThisB1Scale ./ repmat( normSims, [1 size(t2SimSigs4ThisB1Scale,2)] );

      simDataDiff = normdSimSigs - repmat( normdVec', [size(normdSimSigs,1) 1] );
      normDiff = norms( simDataDiff, 2, 2 );
      [~,minDiffIndx] = min( normDiff );
      thisMapCol(j) = T2s( minDiffIndx );   %#ok<PFBNS>
    end
    t2Map{ 1, i } = thisMapCol;
  end
  t2Map = cell2mat( t2Map );
  p.clean;

  if numel( mask ) > 0, t2Map( mask==0 ) = 0; end
end

