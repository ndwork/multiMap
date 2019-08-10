
function [t1Map,m0Map] = mri_mapT1( dataCube, TEs, varargin )
  % [t1Map,m0Map] = mri_mapT1( dataCube, TEs, [, 'alg', alg, 'M0', M0, ...
  %   'MzAt0', MzAt0, 'mask', mask, 'verbose', verbose ] )
  %
  % Fits T1 to an exponential recovery.  If the data is complex, uses the
  % magnitude of the data.
  %
  % Inputs:
  % dataCube - a 3D array of size MxNxK.  Each k index corresponds to an
  %   image of Mz taken at time TEs(k)
  % TEs - a 1D array of size K specifying the spin echo times of each image
  %
  % Optional Inputs:
  % alg - either 'lsqr' or 'linear' (linear fits line to log of data)
  % mask - a 2D array of size MxN.  Only map pixels with nonzero mask values.
  % verbose - scalar; info statements made if verbose is nonzero
  %
  % Outputs:
  % t2Map - a 2D array of size MxN; units are same as TEs
  % m0Map - a 2D array of size MxN; values are proportional to M0
  %
  % Written by Nicholas Dwork - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  p = inputParser;
  p.addParameter( 'MzAt0', [], @isnumeric );
  p.addParameter( 'mask', [], @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'verbose', 0, @(x) isnumeric(x) || islogical(x) );
  p.parse( varargin{:} );
  MzAt0 = p.Results.MzAt0;
  mask = p.Results.mask;
  verbose = p.Results.verbose;

  sData = size( dataCube );
  if numel( mask ) == 0, mask=ones( sData(1:2) ); end
  TEs = TEs(:);

  if numel( MzAt0 ) == 0
    MzAt0 = zeros( sData(1:2) );
  end
  
  [t1Map,m0Map] = mri_mapT1M0( dataCube, TEs, MzAt0, mask, verbose );
end


function [t1Map,m0Map] = mri_mapT1M0( dataCube, TEs, MzAt0, mask, verbose )

  fminconOptions = optimoptions('fmincon','Display','off');
  signalModel = @(sigTEs,thisMzAt0,tmp) tmp(1) .* ( 1 - exp( -sigTEs ./ tmp(2) ) ) + ...
    thisMzAt0 .* exp( -sigTEs ./ tmp(2) );
    % tmp = [ M0, T1 ]

  sData = size( dataCube );
  m0MapCols = cell( 1, sData(2) );
  t1MapCols = cell( 1, sData(2) );

  p = parforProgress( sData(2) );
  parfor i=1:sData(2)
    if verbose ~= 0, p.progress( i, 10 ); end    %#ok<PFBNS>

    m0MapCol = zeros( sData(1), 1 );   %#ok<PFBNS>
    t1MapCol = zeros( sData(1), 1 );

    if max( mask(:,i) ) > 0
      for j=1:sData(1)
        if mask(j,i)==0, continue; end   %#ok<PFBNS>

        %sigTEs = [ 0; TEs; ];
        %thisData = [ 0; abs( squeeze( dataCube(j,i,:) ) ) ];
        sigTEs = TEs;
        thisData = abs( squeeze( dataCube(j,i,:) ) );

        thisMzAt0 = MzAt0( j, i );
        x0 = [1; 1000;];
        xLB = [0; 1;];
        params = fmincon( @(tmp) norm( thisData - signalModel(sigTEs,thisMzAt0,tmp) ), ...
          x0, [], [], [], [], xLB, [], [], fminconOptions );
        m0MapCol(j) = params(1);
        t1MapCol(j) = params(2);

        %figure;  plotnice( thisTEs, thisData );
        %hold all;  plotnice( thisTEs, signalModel(thisTEs,params) );
        %legendnice( 'data', 'model fit' );
      end
    end

    m0MapCols{i} = m0MapCol;
    t1MapCols{i} = t1MapCol;
  end
  p.clean;

  m0Map = cell2mat( m0MapCols );
  t1Map = cell2mat( t1MapCols );
end

