
function run_multiMap
  close all; clear;

  datacase = 3;
  showScale = 3;
  fieldStrength = 1.5;  % Tesla
  verbose = 0;

  [recons,acqTimes,rfSatTime,rf30Time,seTimes,TR,tSat] = loadDatacase( ...
    datacase, 'showScale', showScale, 'verbose', verbose );                                %#ok<ASGLU>

  rfSatDur = 3.2;  % time in ms
  rf1Dur = 3.2;  % tim in ms

  rf180_1Time = seTimes(1);
  rf180_2Time = seTimes(2);

  nSegs = 2;
  nImgTypes = size( recons, 4 ) / nSegs;
  nSlices = size( recons, 3 );


  for sliceIndx=1:nSlices
    sliceRecons = squeeze( recons( :, :, sliceIndx, : ) );
    
    sliceRecons60 = sliceRecons(:,:,1:nImgTypes);
    sliceRecons120 = sliceRecons(:,:,nImgTypes+1:end);

    magMask = max( abs(sliceRecons60(:,:,6)), ...
                   abs(sliceRecons120(:,:,6)) ) > 0.2;
    magMask = imopen( magMask, ones(3) );
    magMask = imclose( magMask, ones(3) );
    magMask = imerode( magMask, ones(5) );
    if verbose ~= 0
      figure('Name','Mask');
      imshowscale( magMask, showScale );
    end


%     %-- B0 Mapping
%     if verbose ~= 0, disp('Working on Off Resonance'); end;
%     offResMap = mri_mapOffResSimple( sliceRecons60(:,:,[1 3]), acqTimes(1:2) ) / (2*pi);
%     figure('Name','Off Resonance (kHz)');
%     imshowscale( magMask .* offResMap, showScale );
%     colorbarnice;  title('Off Resonance (kHz)');
% 
% 
%     %-- Double Angle Method for B1 Mapping
%     b1DataCube = squeeze( sliceRecons( :, :, [6,nImgTypes+6] ) );
%     simple = 1;
%     if simple == 1
%       singleAngleImg = abs( sliceRecons60(:,:,6) );
%       doubleAngleImg = abs( sliceRecons120(:,:,6) );
%       angleMap = acos( doubleAngleImg ./ ( 2 * singleAngleImg ) );
%       b1ScaleMap = angleMap * 180/pi / 60;
%     else
%       b1ScaleMap = mri_mapB1( b1DataCube );
%     end
%     b1Angle60_deg = b1ScaleMap * 60;
%     %if verbose ~= 0
%       figure;  imshowscale( magMask .* b1Angle60_deg, showScale );
%       colorbarnice;  titlenice('B1 Map (degrees)' );  drawnow;
%     %end
% 
% 
%     %-- Correct data with B1 map
%     for imgIndx = 1:4
%       sliceRecons60(:,:,imgIndx) = sliceRecons60(:,:,imgIndx) ./ ...
%         sin( b1ScaleMap * pi/2 ) * sin( pi/2 );
%       sliceRecons120(:,:,imgIndx) = sliceRecons120(:,:,imgIndx) ./ ...
%         sin( b1ScaleMap * pi/2 ) * sin( pi/2 );
%     end
%     sliceRecons60(:,:,5) = sliceRecons60(:,:,5) ./ ...
%       sin( b1ScaleMap * 30*pi/180) * sin( 30*pi/180 );
%     sliceRecons120(:,:,5) = sliceRecons120(:,:,5) ./ ...
%       sin( b1ScaleMap * 30*pi/180 ) * sin( 30*pi/180 );
%     for imgIndx = 6:8
%       sliceRecons60(:,:,imgIndx) = sliceRecons60(:,:,imgIndx) ./ ...
%         sin( b1ScaleMap * 60*pi/180 ) * sin( 60*pi/180 );
%       sliceRecons120(:,:,imgIndx) = sliceRecons120(:,:,imgIndx) ./ ...
%         sin( b1ScaleMap * 120*pi/180 ) * sin( 120*pi/180 );
%     end
% 
% 
% 
%     %-- T1 Mapping
%     if verbose ~= 0, disp('Working on T1 Over M0'); end;
%     I1 = squeeze( sliceRecons60(:,:,5) );  time1 = acqTimes(5);
%     I2 = squeeze( sliceRecons60(:,:,6) );  time2 = acqTimes(6);
%     t1OverM0Img = estimateT1OverM0( I1, I2, time1, time2, ...
%       30*pi/180, 60*pi/180, magMask );
%     
%     %if verbose ~= 0
%       figure;
%       imshowscale( magMask .* t1OverM0Img, showScale, 'range', [0 500] );
%       colorbarnice;  titlenice('M0 Over T1');
%     %end
%
%
%     %-- T2 Star Imaging
%     t2StarData = squeeze( sliceRecons60( :, :, [1:4] ) );
%     TEs = acqTimes([1:4]) - rfSatDur/2;
%     if verbose ~= 0, disp('Working on T2* imaging'); end;
%     t2StarMap = mri_mapT2_linear( t2StarData, TEs, 'mask', magMask, ...
%       'verbose', 1 );
%     t2StarMap( abs(t2StarMap) > 100 ) = 10000;
% %save( 't2StarMap.mat', 't2StarMap' );
% %load t2starMap.mat
%     if verbose ~= 0
%       figure;  imshowscale( magMask .* t2StarMap, showScale, 'range', [0 100] );
%       colorbarnice;  titlenice('T2* Map');  drawnow;
%     end
% 
% 
%     %-- T2 Imaging
%     t2Data = squeeze( sliceRecons60( :, :, 6:8 ) );
%     TEs = acqTimes(6:8) - (tSat + rf1Dur/2);
%     if verbose ~= 0, disp('Working on T2 imaging'); end;
%     %t2Map = mri_mapT2( t2Data, TEs, 'mask', magMask, 'verbose', verbose );
%     [wImg,fImg,t2Map,dB0Map] = mri_wfMultiEchoFit( t2Data, TEs, ...
%       fieldStrength, 'mask', magMask );
% %save( 't2Map.mat', 't2Map' );
% %load t2Map.mat
%     if verbose ~= 0
%       figure; imshowscale( magMask .* t2Map, showScale, 'range', [0 30] );
%       colorbarnice;  titlenice('T2 Map');  drawnow;
%     end


    %-- Water / Fat Fraction with IDEAL
    %[wImg,fImg] = mri_multiEchoFit( sliceRecons60(:,:,1:4), acqTimes(1:4), ...
    %  fieldStrength, 'mask', magMask, 'offResMap', offResMap, 't2Star', t2StarMap );
    %[wImg,fImg,t2StarMap] = mri_multiEchoFit( sliceRecons60(:,:,1:4), acqTimes(1:4), ...
    %  fieldStrength, 'mask', magMask, 'offResMap', offResMap );
    teTimes = acqTimes(1:4) - rfSatTime;
    %[wImg,fImg,t2StarMap,dB0Map] = mri_wfMultiEchoFit_v1( sliceRecons60(:,:,1:4), ...
    %  teTimes, fieldStrength, 'mask', magMask );
    [wMap,fMap,t2StarWMap,t2StarFMap,db0Map] = mri_wfMultiEchoFit( ...
      sliceRecons60(:,:,1:4), teTimes, fieldStrength, 'mask', magMask );
    save( 'ideal.mat', 'wImg', 'fImg', 't2StarMap' );
    % TODO:  fit to 1:4 of sliceRecons60 and sliceRecons120 for added SNR

verbose = 1;
    if verbose ~= 0
      figure;  imshowscale( magMask .* abs(wImg), showScale, 'range', 'nice' );
      colorbarnice;  titlenice('Water Image');

      figure;  imshowscale( magMask .* abs(fImg), showScale, 'range', 'nice' );
      colorbarnice;  titlenice('Fat Image');

      fatFraction = abs( fImg ) ./ ( abs(fImg) + abs(wImg) );
      figure;  imshowscale( magMask .* fatFraction, showScale );
      colorbarnice;  titlenice('Fat Fraction');

      waterFraction = abs( wImg ) ./ ( abs(fImg) + abs(wImg) );
      figure;  imshowscale( magMask .* waterFraction, showScale );
      colorbarnice;  titlenice('Water Fraction');

      figure;  imshowscale( magMask .* t2StarMap, showScale, 'range', [0 50] );
      colorbarnice;  titlenice('T2*');

      figure;  imshowscale( magMask .* dB0Map, showScale );
      colorbarnice;  title('dB0');
    end

    disp('I got here');
  end

end



