function dim_num = hammersley_dim_num_get ( dummy )

%% HAMMERSLEY_DIM_NUM_GET gets the dimension of the leaped Hammersley subsequence.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    18 July 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, dummy, is a dummy argument to be deleted once I
%    figure out how to convince MATLAB to allow me to declare
%    a function with no arguments, but whose invocation can
%    include parentheses!
%
%    Output, integer DIM_NUM, the current value of the dimension.   
%
  global hammersley_DIM_NUM

  dim_num = hammersley_DIM_NUM;

  return
end
