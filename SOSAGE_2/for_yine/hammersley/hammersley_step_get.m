function step = hammersley_step_get ( dummy )

%% HAMMERSLEY_STEP_GET gets the "step" of the leaped Hammersley subsequence.
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
%    Input, DUMMY, is a dummy argument to be deleted once I
%    figure out how to convince MATLAB to allow me to declare
%    a function with no arguments, but whose invocation can
%    include parentheses!
%
%    Output, integer STEP, the step of the leaped Hammersley subsequence.
%
  global hammersley_STEP
%
  step = hammersley_STEP;

  return
end
