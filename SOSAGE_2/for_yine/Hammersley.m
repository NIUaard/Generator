function r = Hammersley(Image,Number)
  dim_num = 4;
  n = Number;
  step = 1;
  nmax=1000;
  seed(1:dim_num) = 0;
  leap(1:dim_num) = 1;
  base(1) = -nmax;
  for i = 2 : dim_num
    base(i) = prime ( i - 1 );
  end

  hammersley_dim_num_set ( dim_num );
  hammersley_step_set ( step );
  hammersley_seed_set ( seed );
  hammersley_leap_set ( leap );
  hammersley_base_set ( base );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  DIM_NUM = %12d\n', dim_num );
  fprintf ( 1, '  N =    %12d\n', n );
  fprintf ( 1, '  STEP = %12d\n', step );
  i4vec_transpose_print ( dim_num, seed, '  SEED = ' );
  i4vec_transpose_print ( dim_num, leap, '  LEAP = ' );
  i4vec_transpose_print ( dim_num, base, '  BASE = ' );

  r = hammersley_sequence ( n );
  for i=1:Number
    M(i)=1;
  end;
  
