
function [M0Img,T1Img] = estimateM0T1Img_simple( I1, I2, ...
  time1, time2, alpha1, alpha2, magMask )

  function out = fModel( x )
    M0=x(1);
    T1=x(2);
    model1 = M0 * ( time1 / T1 ) * sin( alpha1 );
    model2 = M0 * ( 1 - exp(-time2/T1) ) * sin( alpha2 );
    out = abs( model1 - data1 )^2 + abs( model2 - data2 )^2;
  end

  sImg = size( I1 );
  M0Img = zeros( sImg );
  T1Img = zeros( sImg );
  fminConOptions = optimoptions('fmincon','Display','off');
  for i=1:sImg(2)
    if mod(i,10)==0, disp([ 'Working on ', num2str(i), ' of ', num2str(sImg(2)) ]); end;
    for j=1:sImg(1)
      if magMask(j,i) == 0, continue; end;
        x0 = [ 1, 1000 ];
        data1 = abs( I1(j,i) );
        data2 = abs( I2(j,i) );
        bestX = fmincon( @fModel, x0, [], [], [], [], [0, 1d-8], [], [], ...
          fminConOptions );
        M0Img(j,i) = bestX(1);
        T1Img(j,i) = bestX(2);
    end
  end

end

