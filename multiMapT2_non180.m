
function [t2Map,m0Map] = multiMapT2_non180( dataCube1, dataCube2, b1ScaleMap, TEs, varargin )
  % [t2Map,m0Map] = mri_mapT2( dataCube, TEs, [, 'alg', alg, ...
  %   'mask', mask, 'verbose', verbose ] )
  %
  % Fits T2 to an exponential decay.  If the data is complex, uses the
  % magnitude of the data.
  %
  % Inputs:
  % dataCube - a 3D array of size MxNxK.  Each k index corresponds to an
  %   image taken with a different spin echo time.
  % TEs - a 1D array of size K specifying the spin echo times of each image
  %
  % Optional Inputs:
  % alg - either 'lsqr' or 'linear' (linear fits line to log of data)
  % mask - a 2D array of size MxN.  Only map pixels with nonzero mask values.
  % verbose - scalar; info statements made if verbose is nonzero
  %
  % Outputs:
  % t2Map - a 2D array of size MxN; units are same as TEs
  % m0Map - a 2D array of size MxN; values are proportional to M(0)
  %  Note: M(0) is only related to proton density M0 if in equilibrium at time 0
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'mask', [], @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'verbose', 0, @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  mask = p.Results.mask;
  verbose = p.Results.verbose;

  sData = size( dataCube1 );
  if numel( mask ) == 0, mask=ones( sData(1:2) ); end
  TEs = TEs(:);

  % Get an initial guess with linear fit
  [t2MapGuess,m0MapGuess] = mri_mapT2( dataCube1, TEs, 'alg', 'linear', 'mask', mask );
  
  fminconOptions = optimoptions('fmincon','Display','off');
  t2MapCols = cell( sData(2), 1 );
  m0MapCols = cell( sData(2), 1 );
  p = parforProgress( sData(2) );
  parfor i=1:sData(2)
    if verbose ~= 0, p.progress( i, 10 ); end  %#ok<PFBNS>

    m0MapCol = zeros( sData(1), 1 );  %#ok<PFBNS>
    t2MapCol = zeros( sData(1), 1 );
    for j=1:sData(1)
      if mask(j,i)==0, continue; end
      thisData1 = abs( squeeze( dataCube1(j,i,:) ) );
      thisData2 = abs( squeeze( dataCube2(j,i,:) ) );
      spinAngle = b1ScaleMap(j,i) * pi;
      params0 = [m0MapGuess(j,i); t2MapGuess(j,i);] ;
      params = fmincon( ...
        @(tmp) norm( thisData1 - signalModel(TEs,tmp,spinAngle) ) + ...
               norm( thisData2 - signalModel(TEs,tmp,spinAngle) ), ...
          params0, [], [], [], [], [0; 0;], ...
          [], [], fminconOptions );
      m0MapCol(j) = params(1);
      t2MapCol(j) = params(2);
      
      %figure; plotnice( TEs, thisData1 )
      %hold all; plotnice( TEs, signalModel(TEs,params,spinAngle) )
    end

    m0MapCols{i} = m0MapCol;
    t2MapCols{i} = t2MapCol;
  end
  p.clean;

  m0Map = zeros( sData(1:2) );
  t2Map = zeros( sData(1:2) );
  for i=1:sData(2)
    m0Map(:,i) = m0MapCols{i};
    t2Map(:,i) = t2MapCols{i};
  end
  
end


function out = signalModel( TEs, params, angle )
  % This objective function is taken from "Errors in the Measurements of T2 Using
  % Multiple-Echo MRI Techniques by Majumdar et al. (1986)
  M0 = params(1);
  T2 = params(2);

  ca = cos(angle);
  sa = sin(angle);
  saSq = sa * sa;
  tmp1 = 1 - ca;
  tmp1Sq = tmp1 * tmp1;
  tmp2 = 1 + ca;

  nTEs = numel( TEs );
  out = zeros( numel(TEs), 1 );
  out(1) = 0.5 * tmp1;
  out(2) = 0.25 * tmp1*tmp1;
  if nTEs > 2
    out(3) = 0.125 * ( tmp1*tmp1Sq + tmp1*tmp2*tmp2 ) + 0.5 * ca*saSq;
  end
  if nTEs > 3
    out(4) = 0.0625 * tmp1Sq*tmp1Sq + 0.125 * 3 * saSq*saSq + 0.5*ca*sa*sa*tmp1;
  end
  if nTEs > 4, error( 'Not yet implementes' ); end

  out = out .* M0 .* exp( -TEs ./ T2 );
end
