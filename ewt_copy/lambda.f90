implicit none
integer :: i, j
double reecision :: t

open(unit = 100, file = 'readfile/fort.50', status = 'old', action = 'read')
