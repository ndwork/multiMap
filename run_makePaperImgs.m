
function run_makePaperImgs
  close all; clear; rng(1);

  datacases = 1:16;
  verbose = 0;
  outDir = './outputs/';

  mkdir( outDir );


  run_standardScans;

  for datacase = datacases
    thisIndx = indx2str( datacase, max(datacases) );
    disp([ 'Working on ', thisIndx, ' of ', num2str( max(datacases) ) ]);

    subOutDir = [ '/output_', thisIndx ];

    run_multiMap( datacase, 'outDir', [outDir,subOutDir], 'verbose', verbose );
  end

end
