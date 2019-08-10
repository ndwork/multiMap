
function [db0Map, b1ScaleMap, m0Map, t1Map, t2Map, t2StarMapW, t2StarMapF, ffMap] = multiMap( ...
  recons, acqTimes, rfSatTime, tSat, sliceThickness, seTimes, magThresh, fieldStrength, varargin )
  % [db0Map, b1ScaleMap, m0Map, t1Map, t1OverM0Map, t2Map, t2StarMap, ffMap] = multiMap( ...
  %   recons, acqTimes, rfSatTime, tSat, sliceThickness, seTimes, magThresh, fieldStrength [, ...
  %   'b0Bound', b0Bound, 'showScale', showScale, 'verbose', true/false )
  %
  % Inputs:
  % sliceThickness - slice thickness used for RF pulse (in cm)
  % Written by Nicholas Dwork - Copyright 2019
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'b0Bound', 0.8, @ispositive );
  p.addParameter( 'showScale', 5, @ispositive );
  p.addParameter( 'verbose', false, @(x) islogical(x) || isnumeric(x) );
  p.parse( varargin{:} );
  b0Bound = p.Results.b0Bound;
  showScale = p.Results.showScale;
  verbose = p.Results.verbose;
  excFlipAngle = 60;
  fracThresh = 0.05;
  includeMzAt0 = true;

  nSegs = 2;
  nImgTypes = size( recons, 4 ) / nSegs;
  nSlices = size( recons, 3 );

  if verbose == true
    figure; imshowscale( squeeze( abs( recons( :, :, 1, 5 ) ) ), showScale );
    titlenice('Before cleaning');  drawnow;
  end
  %recons = cleanRecons( recons );
  %if verbose == true
  %  figure; imshowscale( squeeze( abs( recons( :, :, 1, 5 ) ) ), showScale );
  %  titlenice( 'After cleaning' );
  %end


  for sliceIndx=1:nSlices
    sliceRecons = squeeze( recons( :, :, sliceIndx, : ) );

    mask = makeMaskFromRecons( sliceRecons, magThresh );
    if verbose == true
      figure; imshowscale( mask, showScale );  titlenice( 'Mask' );
    end

    sliceRecons60 = sliceRecons(:,:,1:nImgTypes);
    sliceRecons120 = sliceRecons(:,:,nImgTypes+1:end);


    %-- Off Resonance Mapping
    if verbose ~= 0, disp('Working on Off Resonance'); end
    offResIndx2nd = 3;
    offResMap_kHz = mapOffRes( sliceRecons60(:,:,[1 offResIndx2nd]), ...
      acqTimes([1,offResIndx2nd]), 'mask', mask ) / (2*pi);
    gammaBar = getGammaH;
    db0Map = offResMap_kHz / gammaBar / 10000;
    if verbose == true
      figure;  imshowscale( mask .* db0Map, showScale, 'range', 'nice' );
      colorbarnice;  titlenice('Off Resonance (kHz)');  drawnow;
    end


    %-- Double Angle Method for B1 Mapping
    if verbose ~= 0, disp('Working on B1'); end
    b1DataIndx = nImgTypes-3;
    b1DataCube = squeeze( sliceRecons( :, :, [b1DataIndx,nImgTypes+b1DataIndx] ) );
    b1ScaleMap = mri_doubleAngleMapB1( b1DataCube, sliceThickness, 'mask', mask );
    if verbose == true
      figure;  imshowscale( mask .* b1ScaleMap, showScale, 'range', [0.6 1.3] );
      colorbarnice;  titlenice('B1 Scale Map' );  drawnow;
    end


    %-- T2 Imaging
    if verbose ~= 0, disp('Working on T2'); end
    t2Data60  = squeeze(  sliceRecons60( :, :, nImgTypes-2:end ) );
    t2Data120 = squeeze( sliceRecons120( :, :, nImgTypes-2:end ) );
    TEs = acqTimes(nImgTypes-2:end) - tSat;
    if verbose ~= 0, disp('Working on T2 imaging'); end
    alg = 'lsqr';
    if strcmp( alg, 'linear' )  % linear fit
      t2Map = mri_mapT2( t2Data60, TEs, 'mask', mask, 'alg', 'linear', 'verbose', verbose );
    elseif strcmp( alg, 'lsqrExp' )  % least-squares to exponential
      t2Map = multiMapT2( t2Data60, t2Data120, b1ScaleMap, TEs, ...
        'mask', mask, 'verbose', verbose );
    elseif strcmp( alg, 'lsqr' )  % least-squares to exponential
      t2Map= multiMapT2_non180( t2Data60, t2Data120, b1ScaleMap, TEs, ...
        'mask', mask, 'verbose', verbose );
    elseif strcmp( alg, 'fp' )  % fingerprinting
      t2SpinTimes = seTimes - tSat;
      t2Map = multiMapT2_fp( t2Data60, TEs, t2SpinTimes, excFlipAngle, b1ScaleMap, ...
        'mask', mask, 'verbose', verbose );
    end
    if verbose == true
      figure;  showImageCube( abs(t2Data60), showScale, 'border', 4, 'borderValue', 'max' );
      titlenice( 't2Data60' );  drawnow;
      %figure;  showImageCube( abs(t2Data120), showScale, 'border', 4, 'borderValue', 'max' );
      %titlenice('t2Data120');
      figure; imshowscale( mask .* t2Map, showScale, 'range', [0 50] );
      colorbarnice;  titlenice('T2 Map');  drawnow;  drawnow;
    end


    %-- Water / Fat Fraction Mapping with IDEAL
    nImgs4MultiEchoFit = 5;
    teTimes = acqTimes(1:nImgs4MultiEchoFit) - rfSatTime;
    %[wMap,fMap,t2StarMap,db0Map,df0Map] = mri_wfMultiEchoFit( ...
    %  sliceRecons60(:,:,1:nImgs4MultiEchoFit), teTimes, fieldStrength, ...
    %  'mask', mask, 'offResMap_kHz', offResMap_kHz, 'b0Bound', b0Bound );  %#ok<ASGLU>
    [wMap,fMap,t2StarMapW,t2StarMapF,dB0Map] = mri_wfMultiEchoFit_2T2s( ...
      sliceRecons60(:,:,1:nImgs4MultiEchoFit), teTimes, fieldStrength, ...
      'mask', mask, 'offResMap_kHz', offResMap_kHz, 'b0Bound', b0Bound );  %#ok<ASGLU>
    ffMap = abs(fMap) ./ ( abs(wMap) + abs(fMap) );  % Fat Fraction Map
    ffMap( mask == 0 ) = 0;
    wfMap = abs(wMap) ./ ( abs(wMap) + abs(fMap) );  % Water Fraction Map
    wfMap( mask == 0 ) = 0;
    if verbose == true
      figure; imshowscale( ffMap, showScale, 'range', [0 1.0] );
      colorbarnice;  titlenice( 'Fat Fraction' );  drawnow;
      figure; imshowscale( wfMap, showScale, 'range', [0 1.0] );
      colorbarnice;  titlenice( 'Water Fraction' );  drawnow;
      if exist( 't2StarMap', 'var' )
        figure; imshowscale( t2StarMap, showScale, 'range', [0 50] );
        colorbarnice;  titlenice( 'T2* Map' );  drawnow;
        t2StarMapW = t2StarMap;
        t2StarMapF = t2StarMap;
      elseif exist( 't2StarMapW', 'var' )
        figure; imshowscale( t2StarMapW .* ( wfMap > fracThresh ) , showScale, 'range', [0 50] );
        colorbarnice;  titlenice( 'Water T2* Map' );  drawnow;
        figure; imshowscale( t2StarMapF .* ( ffMap > fracThresh ), showScale, 'range', [0 50] );
        colorbarnice;  titlenice( 'Fat T2* Map' );  drawnow;
      end
    end


    %-- T1 Mapping
    t1Times = [ acqTimes( nImgs4MultiEchoFit + 1 : nImgs4MultiEchoFit + 3 ) ];
    t1Thetas = [ 30, 30, 60 ] * pi/180;
    t1Data60 = squeeze( abs( squeeze( sliceRecons60(:,:,...
      nImgs4MultiEchoFit+1:nImgs4MultiEchoFit+3 ) ) ) );
    t1Data120 = squeeze( abs( squeeze( sliceRecons120(:,:,...
      nImgs4MultiEchoFit+1:nImgs4MultiEchoFit+3 ) ) ) );
    MzAt0 = [];
    if includeMzAt0 == true
      %I0 = abs( squeeze( sliceRecons60(:,:,1) ) );
      I0 = abs( wMap ) + abs( fMap );
      MzAt0 = abs(I0) ./ tan( b1ScaleMap .* pi/2 );
      MzAt0( mask==0 | b1ScaleMap == 0 ) = 0;
      %t1Data = cat( 3, MzAt0, t1Data );
      %t1Thetas = [ pi/2 t1Thetas ];   %#ok<AGROW>
      %t1Times = [ 0 t1Times ];   %#ok<AGROW>
    end
    [t1Map,m0Map] = mri_mapT1_repeatedExcitations( t1Data60, t1Times, t1Thetas, ...
      sliceThickness, 'b1ScaleMap', b1ScaleMap, 'MzAt0', MzAt0, 'verbose', verbose );

    if verbose == true
     figure;  imshowscale( mask .* t1Map, showScale, 'range', [0 300] );
     colorbarnice;  titlenice('T1');  drawnow;
     figure;  imshowscale( mask .* m0Map, showScale, 'range', [0 10] );
     colorbarnice;  titlenice('M0');  drawnow;
    end

    %if verbose ~= 0, disp('Working on M0 Over T1'); end
    %I1 = squeeze( sliceRecons60(:,:,4) );  time1 = acqTimes(4);
    %I2 = squeeze( sliceRecons60(:,:,5) );  time2 = acqTimes(5);
    %m0R1Img = estimateM0R1( I1, I2, time1, time2, ...
    %  30*pi/180, 60*pi/180, 'magMask', mask, 'b1ScaleMap', b1ScaleMap );
    %m0R1Img( mask == 0 ) = 0;
    %if verbose == true
    %  figure;  imshowscale( m0R1Img, showScale, 'range', 'nice' );
    %  colorbarnice;  titlenice('M0 R1');  drawnow;
    %end

    % %-- T1/M0 Mapping
    %t1OverM0Map = estimateT1OverM0( I1, I2, time1, time2, ...
    %  30*pi/180, 60*pi/180, 'magMask', mask, 'b1ScaleMap', b1ScaleMap );
    %if verbose == true
    %  figure;  imshowscale( mask .* t1OverM0Map, showScale, 'range', [0 1000] );
    %  colorbarnice;  titlenice('T1 Over M0');
    %end

  end

end
