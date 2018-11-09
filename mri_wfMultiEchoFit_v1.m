
function [wMap,fMap,t2StarMap,db0Map] = mri_wfMultiEchoFit_v1( recons, ...
  teTimes, fieldStrength, varargin )
  % [wMap,fMap,b0Map] = mri_wfMultiEchoFit( recons, acqTimes, fieldStrength )
  %
  % Based on the 2004 paper "Multicoil Dixon Chemical Species Separation
  % With an Iterative Least-Squares Estimation Method" by Reeder et al.
  % (including Gold and Pelc)
  %
  % Inputs:
  % recons - an M x N x Nacq array
  %   Each reconstructed image is M x N
  %   Nacq is the number of acquisitions to use in the minimization
  % acqTimes - the acquisition times of each recon (in ms)
  % fieldStrength - value of B0 in Tesla
  %
  % Optional Inputs:
  %
  % Outputs:
  % wMap - the water image
  % fMap - the fat image
  % t2StarMap - the t2Star value of each pixel
  % b0Map - the off resonance image
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

%showImageCube( abs( recons ), 3, 'nImgsPerRow', 4 );  titlenice( 'Magnitude' );
%showImageCube( angle( recons ), 3, 'nImgsPerRow', 4 );  titlenice( 'Phase' );

  p = inputParser;
  p.addParameter( 'mask', [] );
  p.parse( varargin{:} );
  mask = p.Results.mask;

  sRecons = size( recons );
  if numel( mask ) == 0, mask = ones( sRecons(1:2) ); end;

  wMapCells = cell( 1, sRecons(2) );
  fMapCells = cell( 1, sRecons(2) );
  db0MapCells = cell( 1, sRecons(2) );
  t2StarMapCells = cell( 1, sRecons(2) );

  p = parforProgress( sRecons(2) );
  parfor n=1:sRecons(2)
%for n=179  % fat
%for n=140  % muscle
%for n=70;  % fat
%for n=69;  % fat in knee corner
%for n = 33;  % left bottle
    p.progress( n, 1 );                                                                    %#ok<PFBNS>

    theseW = zeros( sRecons(1), 1 );                                                       %#ok<PFBNS>
    theseF = zeros( sRecons(1), 1 );
    theseDB0 = zeros( sRecons(1), 1 );
    theseT2Star = zeros( sRecons(1), 1 );
    for m=1:sRecons(1)
%for m=176  % fat
%for m=198  % muscle
%for m=189;  % fat
%for m=101;  % fat in knee corner
%for m=87;  % left bottle
%figure; showFeaturesOnImg( [n,m], abs(recons(:,:,1)), 'scale', 4 );
      if mask(m,n)==0, continue; end;

% w = 0.2;  f=0.8;  df0=0.08;  t2Star=20;
% disp([ w f ]);
% TEs = teTimes;
% gammaBar = getGammaH;  % kHz / Gauss
% fatFreq = fieldStrength * 1d4 * -gammaBar * 3.4d-6;  % kHz
% simSignal = exp( 1i*2*pi*df0*TEs ) .* exp(-TEs/t2Star) .* ( ...
%   w + f * exp( 1i*2*pi*fatFreq*TEs ) )  ;
% [w,f,t2Star,dB0] = multiEchoFit( simSignal, teTimes, fieldStrength );

      thisRecon = squeeze( recons(m,n,:) );
      [w,f,t2Star,dB0] = multiEchoFit( thisRecon, teTimes, fieldStrength );
      theseW(m) = w;   theseDB0(m) = dB0;
      theseF(m) = f;   theseT2Star(m) = t2Star;
    end
    wMapCells{n} = theseW;
    fMapCells{n} = theseF;
    db0MapCells{n} = theseDB0;
    t2StarMapCells{n} = theseT2Star;
  end
  p.clean;
  wMap = cell2mat( wMapCells );  db0Map = cell2mat( db0MapCells );
  fMap = cell2mat( fMapCells );  t2StarMap = cell2mat( t2StarMapCells );
end


function [w,f,t2Star,dB0] = multiEchoFit( data, TEs, fieldStrength )
  gammaBar = getGammaH;  % kHz / Gauss
  fatFreq = fieldStrength * 1d4 * -gammaBar * 3.4d-6;  % kHz
    % Chemical shift of fat is 3.4 ppm.  1T = 10^4 Gauss.
    % Fat precesses slower than water

  dB0s = linspace( -fieldStrength, fieldStrength, 10001 ) * 1.2d-6;
  t2Stars = 1 : 1 : 40;

  df0s = dB0s * 1d4 * -gammaBar;  % 1T = 10^4 Gauss; kHz
  minCost = Inf;
  eFat = exp( 1i * 2*pi * fatFreq * TEs(:) );
  allCosts = zeros( numel( t2Stars ), numel( df0s ) );
  unOffRes = exp( -1i * 2*pi * TEs(:) * df0s(:)' );
  rotSig = bsxfun( @times, data(:), unOffRes );
  for t2StarIndx = 1 : numel( t2Stars )
%for t2StarIndx = 20
    e2Stars = exp( -TEs(:) / t2Stars( t2StarIndx ) );
    A = [ e2Stars, eFat .* e2Stars ];
    rhos = A \ rotSig;
    diffs = A*rhos - rotSig;
    costs = sum( diffs .* conj(diffs), 1 );
    allCosts( t2StarIndx, : ) = costs;
    [thisMinCost,thisMinCostIndx] = min( costs );
    if thisMinCost < minCost
      minCost = thisMinCost;
      minCosts = costs;
      bestDB0 = dB0s( thisMinCostIndx );
      bestDf0 = df0s( thisMinCostIndx );
      bestRho = rhos( :, thisMinCostIndx );
      bestT2Star = t2Stars( t2StarIndx );
    end
  end

  t2Star = bestT2Star;
  dB0 = bestDB0;
  df0 = dB0 * 1d4 * -gammaBar;  % 1T = 10^4 Gauss; kHz
  w = bestRho(1);
  f = bestRho(2);

%   simSignal = exp( 1i*2*pi*df0*TEs ) .* exp(-TEs/t2Star) .* ( ...
%    w + f * exp( 1i*2*pi*fatFreq*TEs ) )  ;
%   figure; stemnice( abs( data ), 'b' ); hold on; plotnice( abs( simSignal ), 'r' );
%   titlenice( 'Magnitude' );
%   figure; stemnice( unwrap( angle( data ) ), 'b' );
%   hold on; plotnice( unwrap( angle( simSignal ) ), 'r' );
%   titlenice( 'Phase' );
%   figure; stemnice( real( data ), 'b' ); hold on; plotnice( real( simSignal ), 'r' );
%   titlenice( 'Real' );
%   figure; stemnice( imag( data ), 'b' ); hold on; plotnice( imag( simSignal ), 'r' );
%   titlenice( 'Imag' );
%   figure; plotnice( df0s, minCosts );  titlenice( 'minCosts' );
%   figure; imshowscale( allCosts, 3 ); titlenice( 'allCosts' );
end
