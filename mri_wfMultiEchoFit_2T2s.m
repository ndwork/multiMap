
function [wMap,fMap,t2StarMapW,t2StarMapF,db0Map,df0Map] = mri_wfMultiEchoFit_2T2s( ...
  recons, teTimes, fieldStrength, varargin )
  % [wMap,fMap,t2StarMapW,t2StarMapF,db0Map,df0Map] = mri_wfMultiEchoFit( ...
  %   recons, teTimes, fieldStrength )
  %
  % Based on the 2004 paper "Multicoil Dixon Chemical Species Separation
  % With an Iterative Least-Squares Estimation Method" by Reeder et al.
  % (including Gold and Pelc)
  %
  % Inputs:
  % recons - an M x N x Nacq array
  %   Each reconstructed image is M x N
  %   Nacq is the number of acquisitions to use in the minimization
  % teTimes - the echo times of each recon (in ms)
  % fieldStrength - value of B0 in Tesla
  %
  % Optional Inputs:
  % b0Bound - in microTesla (default is 0.5)
  % mask - sets all values where mask is 0 to 0
  % offResMap_kHz - off-resonance map; provides an initial guess for db0
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
  p.addParameter( 'b0Bound', 0.5, @isnumeric );
  p.addParameter( 'mask', [] );
  p.addParameter( 'offResMap_kHz', [], @isnumeric );
  p.parse( varargin{:} );
  b0Bound = p.Results.b0Bound;
  mask = p.Results.mask;
  offResMap_kHz = p.Results.offResMap_kHz;

  sRecons = size( recons );
  if numel( mask ) == 0, mask = ones( sRecons(1:2) ); end

  if numel( offResMap_kHz ) > 0
    db0_guess = offResMap_kHz * 1d-4 / -getGammaH();
  else
    db0_guess = zeros( sRecons(1:2) );
  end

  wMapCells = cell( 1, sRecons(2) );
  fMapCells = cell( 1, sRecons(2) );
  db0MapCells = cell( 1, sRecons(2) );
  df0MapCells = cell( 1, sRecons(2) );
  t2StarMapCellsW = cell( 1, sRecons(2) );
  t2StarMapCellsF = cell( 1, sRecons(2) );

  p = parforProgress( sRecons(2) );
  parfor n=1:sRecons(2)
    p.progress( n, 1 );   %#ok<PFBNS>

    theseW = zeros( sRecons(1), 1 );   %#ok<PFBNS>
    theseF = zeros( sRecons(1), 1 );
    theseDB0 = zeros( sRecons(1), 1 );
    theseDf0 = zeros( sRecons(1), 1 );
    theseT2StarW = zeros( sRecons(1), 1 );
    theseT2StarF = zeros( sRecons(1), 1 );

    for m=1:sRecons(1)
      if mask(m,n)==0, continue; end

      thisRecon = squeeze( recons(m,n,:) );
      [w,f,t2StarW,t2StarF,dB0,df0] = multiEchoFit( thisRecon(:), teTimes(:), ...  %#ok<PFBNS>
        fieldStrength, db0_guess(m,n), b0Bound );
      theseW(m) = w;   theseF(m) = f;
      theseDB0(m) = dB0;  theseDf0(m) = df0;
      theseT2StarW(m) = t2StarW;  theseT2StarF(m) = t2StarF;
    end
    wMapCells{n} = theseW;
    fMapCells{n} = theseF;
    db0MapCells{n} = theseDB0;
    df0MapCells{n} = theseDf0;
    t2StarMapCellsW{n} = theseT2StarW;
    t2StarMapCellsF{n} = theseT2StarF;
  end
  p.clean;
  wMap = cell2mat( wMapCells );
  db0Map = cell2mat( db0MapCells );
  df0Map = cell2mat( df0MapCells );
  fMap = cell2mat( fMapCells );
  t2StarMapW = cell2mat( t2StarMapCellsW );
  t2StarMapF = cell2mat( t2StarMapCellsF );
end


function [w,f,t2StarW,t2StarF,dB0,df0] = multiEchoFit( data, TEs, fieldStrength, db0_guess, b0Bound )
  gammaBar = getGammaH;  % kHz / Gauss
  fatFreq = fieldStrength * 1d4 * -gammaBar * 3.4d-6;  % kHz
    % Chemical shift of fat is 3.4 ppm.  1T = 10^4 Gauss.
    % Fat precesses slower than water

  dB0s = linspace( -fieldStrength, fieldStrength, 101 ) * b0Bound * 1d-6;
  dB0s = dB0s + db0_guess;
  t2Stars = 1 : 1 : 100;
  %df0_guess = db0_guess * 1d4 * -gammaBar

  df0s = dB0s * 1d4 * -gammaBar;  % 1T = 10^4 Gauss; kHz
  minCost = Inf;
  eFat = exp( 1i * 2*pi * fatFreq * TEs );
  %allCosts = zeros( numel( t2Stars ), numel( df0s ) );
  unOffRes = exp( -1i * 2*pi * TEs * df0s(:)' );
  rotSig = bsxfun( @times, data, unOffRes );
  for indxT2StarW = 5 : numel( t2Stars )
    e2StarsW = exp( -TEs / t2Stars( indxT2StarW ) );

    for indxT2StarF = 5 : numel( t2Stars )      
      e2StarsF = exp( -TEs / t2Stars( indxT2StarF ) );

      A = [ e2StarsW, eFat .* e2StarsF ];
      rhos = A \ rotSig;
      diffs = A*rhos - rotSig;
      costs = sum( diffs .* conj(diffs), 1 );
      [thisMinCost,thisMinCostIndx] = min( costs );
      if thisMinCost < minCost
        minCost = thisMinCost;
        bestDB0 = dB0s( thisMinCostIndx );
        bestDf0 = df0s( thisMinCostIndx );
        bestRho = rhos( :, thisMinCostIndx );
        bestT2StarW = t2Stars( indxT2StarW );
        bestT2StarF = t2Stars( indxT2StarF );
      end
    end
  end

  t2StarW = bestT2StarW;
  t2StarF = bestT2StarF;
  dB0 = bestDB0;
  df0 = bestDf0;  % 1T = 10^4 Gauss; kHz
  w = bestRho(1);
  f = bestRho(2);

  %bestE2StarsW = exp( -TEs / bestT2StarW );
  %bestE2StarsF = exp( -TEs / bestT2StarF );
  %bestA = [ bestE2StarsW, eFat .* bestE2StarsF ];
  %model = bestA * bestRho;
  %figure; plotnice( TEs, abs(model) );
  %bestRotSig = exp( -1i * 2*pi * TEs * bestDf0 ) .* data;
  %hold on; plotnice( TEs, abs(bestRotSig) );
  %legend( 'model', 'data' );
end

