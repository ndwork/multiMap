
function [wMap,fMap,t2StarWMap,t2StarFMap,db0Map] = mri_wfMultiEchoFit( ...
  recons, teTimes, fieldStrength, varargin )
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
  t2StarWCells = cell( 1, sRecons(2) );
  t2StarFCells = cell( 1, sRecons(2) );

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
    theseT2StarW = zeros( sRecons(1), 1 );
    theseT2StarF = zeros( sRecons(1), 1 );
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
      [w,f,t2StarW,t2StarF,dB0] = multiEchoFit( thisRecon, teTimes, ...
        fieldStrength );
      theseW(m) = w;   theseDB0(m) = dB0;
      theseF(m) = f;   theseT2StarW(m) = t2StarW;  theseT2StarF = t2StarF;
    end
    wMapCells{n} = theseW;
    fMapCells{n} = theseF;
    db0MapCells{n} = theseDB0;
    t2StarWCells{n} = theseT2StarW;
    t2StarFCells{n} = theseT2StarF;
  end
  p.clean;
  wMap = cell2mat( wMapCells );  t2StarWMap = cell2mat( t2StarWCells );
  fMap = cell2mat( fMapCells );  t2StarFMap = cell2mat( t2StarFCells );
  db0Map = cell2mat( db0MapCells );
end


function [w,f,t2StarW,t2StarF,dB0] = multiEchoFit( data, TEs, fieldStrength )
  gammaBar = getGammaH;  % kHz / Gauss
  fatFreq = fieldStrength * 1d4 * -gammaBar * 3.4d-6;  % kHz
    % Chemical shift of fat is 3.4 ppm.  1T = 10^4 Gauss.
    % Fat precesses slower than water

  dB0s = linspace( -fieldStrength, fieldStrength, 10001 ) * 1.2d-6;
  t2StarsW = 1 : 1 : 40;
  t2StarsF = 1 : 1 : 40;

  df0s = dB0s * 1d4 * -gammaBar;  % 1T = 10^4 Gauss; kHz
  minCost = Inf;
  eFat = exp( 1i * 2*pi * fatFreq * TEs(:) );
  unOffRes = exp( -1i * 2*pi * TEs(:) * df0s(:)' );
  rotSig = bsxfun( @times, data(:), unOffRes );
  nTEs = numel( TEs );
  for t2StarWIndx = 1 : numel( t2StarsW )
    e2StarsW = exp( -TEs(:) / t2StarsW( t2StarWIndx ) );

    A = [ e2StarsW, zeros(nTEs,1) ];
    for t2StarFIndx = 1 : numel( t2StarsF )
      e2StarsF = exp( -TEs(:) / t2StarsF( t2StarFIndx ) );
      A(:,2) = eFat .* e2StarsF;
      rhos = A \ rotSig;
      diffs = A*rhos - rotSig;
      costs = sum( diffs .* conj(diffs), 1 );
      [thisMinCost,thisMinCostIndx] = min( costs );
      if thisMinCost < minCost
        minCost = thisMinCost;
        minCosts = costs;
        bestDB0 = dB0s( thisMinCostIndx );
        bestDf0 = df0s( thisMinCostIndx );
        bestRho = rhos( :, thisMinCostIndx );
        bestT2StarW = t2StarsW( t2StarWIndx );
        bestT2StarF = t2StarsF( t2StarFIndx );
      end
    end
  end

  t2StarW = bestT2StarW;
  t2StarF = bestT2StarF;
  dB0 = bestDB0;
  df0 = dB0 * 1d4 * -gammaBar;  % 1T = 10^4 Gauss; kHz
  w = bestRho(1);
  f = bestRho(2);

  simSignal = exp( 1i*2*pi*df0*TEs ) .* ( w .* exp(-TEs/t2StarW) + ...
    f .* exp(-TEs/t2StarW) .* exp( 1i*2*pi*fatFreq*TEs ) )  ;
  figure; stemnice( abs( data ), 'b' ); hold on; plotnice( abs( simSignal ), 'r' );
  titlenice( 'Magnitude' );
  figure; stemnice( unwrap( angle( data ) ), 'b' );
  hold on; plotnice( unwrap( angle( simSignal ) ), 'r' );
  titlenice( 'Phase' );
  figure; stemnice( real( data ), 'b' ); hold on; plotnice( real( simSignal ), 'r' );
  titlenice( 'Real' );
  figure; stemnice( imag( data ), 'b' ); hold on; plotnice( imag( simSignal ), 'r' );
  titlenice( 'Imag' );
  figure; plotnice( df0s, minCosts );  titlenice( 'minCosts' );
end
