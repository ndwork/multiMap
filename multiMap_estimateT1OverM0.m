
function [t1OverM0s,relErrs] = multiMap_estimateT1OverM0( sigTimes, sigThetas, data, varargin )
  % [t1OverM0s,relErrs] = estimateT1OverM0_compound( t1Times, t1Thetas, t1Data [, ...
  %   'mask', mask, 'b1ScaleMap', b1ScaleMap ] )
  %
  % This algorithm first tries to fit linear recovery.  Only after doing so, does it fit
  % an exponential curve

  p = inputParser;
  p.addParameter( 'mask', [], @(x) isnumeric(x) || islogical(x) );
  p.addParameter( 'b1ScaleMap', [], @isnumeric );
  p.parse( varargin{:} );
  mask = p.Results.mask;
  b1ScaleMap = p.Results.b1ScaleMap;

  if numel( b1ScaleMap ) == 0, b1ScaleMap = ones( size(I1) ); end

  [Ny,Nx,nSigs] = size( data );
  b1Scalings = repmat( b1ScaleMap, [ 1 1 nSigs ] );

  thetasForBsxfun = reshape( sigThetas, [ 1 1 nSigs ] );
  scaledThetas = bsxfun( @times, b1Scalings, thetasForBsxfun );
  sinScaledThetas = sin( scaledThetas );

  t1OverM0s = zeros( Ny, Nx );
  sigTimes = sigTimes(:);

  fminconOptions = optimoptions( 'fmincon', 'Display', 'off' );




  % First fit a linear recovery to the data
  relErrCols = cell(1,Nx);
  t1OverM0Cols = cell(1,Nx);
  parfor i = 1 : Nx
%for i=32  % (92,32) fat phantom
%for i=56  % (60,56) middle phantom
    disp([ 'Working on ', num2str(i), ' of ', num2str(Nx) ]);

    relErrCol = zeros(Ny,1);
    t1OverM0Col = zeros(Ny,1);

    for j = 1 : Ny
%for j=92  % (92,32) fat phantom
%for j=60  % (60,56) middle phantom
      if mask(j,i) == 0, continue; end

      theseSinScaledThetas = sinScaledThetas( j, i, : );
      theseSinScaledThetas = theseSinScaledThetas(:);
      theseData = abs( squeeze( data(j,i,:) ) );
      params0 = [0; 0;];
      params = fmincon( @(tmp) norm( ...
        theseData - linearDataModel(tmp, sigTimes, theseSinScaledThetas), 2 ), ...
        params0, [], [], [], [], [-Inf,0], [], [], fminconOptions );
      modelResult = linearDataModel( params, sigTimes, theseSinScaledThetas );
      %figure; plotnice( sigTimes, theseData ); hold on; plotnice( sigTimes, modelResult, 'r' );

      t1OverM0Col(j) = params(2);
      relErrCol(j) = norm( theseData - modelResult, 2 ) / norm( theseData, 2 );

    end
    
    relErrCols{i} = relErrCol;
    t1OverM0Cols{i} = t1OverM0Col;
  end
  t1OverM0s = cell2mat( t1OverM0Cols );
  relErrs = cell2mat( relErrCols );

end


function out = linearDataModel( params, sigTimes, sinScaledThetas )
  modelMzs = params(1) + params(2) * sigTimes;
  out = modelMzs .* sinScaledThetas;
end



