function hammersley_test005 ( )

%% TEST005 compares "classical" Hammersley sequences of different length.
%
%  Discussion:
%
%    The point made by this example is simple.  In a "classical" Hammersley
%    sequence, the first dimension is treated specially.  A sequence length
%    NMAX is specified, and the elements of the first dimension are
%    simply 0/NMAX, 1/NMAX, 2/NMAX, ..., NMAX/NMAX.
%
%    In the second and higher dimensions, a different scheme is used.
%    A classical Hammersley sequence is based on primes, so dimension 2
%    would be associated with 2, dimension 3 with 3, dimension 4 with 5
%    and so on.
%
%    If two classical Hammersley sequences are generated, differing only
%    in the value of NMAX, then one sequence will repeat all the values of
%    the other, except in the first dimension.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    05 May 2008
%
%  Author:
%
%    John Burkardt
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST005\n' );
  fprintf ( 1, '  I4_TO_HAMMERSLEY_SEQUENCE computes N elements of\n' );
  fprintf ( 1, '  a Hammersley sequence on a single call.\n' );
  fprintf ( 1, '  All arguments are specified explicitly.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Compare "classical" Hammersley sequences of length 11\n' );
  fprintf ( 1, '  and 16.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Note that the second sequence repeats the values of the\n' );
  fprintf ( 1, '  first sequence, except in the first dimension.\n' );

  for i = 1 : 3

    if ( i == 1 )
      nmax = 5;
    elseif ( i == 2 )
      nmax = 10;
    elseif ( i == 3 )
      nmax = 15;
    end 

    dim_num = 4;
    n = nmax + 1;
    step = 0;
    seed(1:dim_num) = 0;
    leap(1:dim_num) = 1;
    base(1) = -(nmax);
    for i = 2 : dim_num
      base(i) = prime ( i - 1 );
    end

    fprintf ( 1, '\n' );
    fprintf ( 1, '  DIM_NUM = %12d\n', dim_num );
    fprintf ( 1, '  N =    %12d\n', n );
    fprintf ( 1, '  STEP = %12d\n', step );
    i4vec_transpose_print ( dim_num, seed, '  SEED = ' );
    i4vec_transpose_print ( dim_num, leap, '  LEAP = ' );
    i4vec_transpose_print ( dim_num, base, '  BASE = ' );

    r = i4_to_hammersley_sequence ( dim_num, n, step, seed, leap, base );

    fprintf ( 1, '\n' );
    fprintf ( 1, '    STEP   Hammersley\n' );
    fprintf ( 1, '\n' );
    for j = 1 : n
      fprintf ( 1, '  %6d  ', step+j-1 );
      for i = 1 : dim_num
        fprintf ( 1, '%12f  ', r(i,j) );
      end
      fprintf ( 1, '\n' );
    end

  end

  return
end
