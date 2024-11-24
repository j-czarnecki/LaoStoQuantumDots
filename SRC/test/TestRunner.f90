
! TestRunner.f90 - runs fUnit test suites
!
! funit generated this file on 2024-11-24 15:35:18 +0100.

program TestRunner

    use combinatory_fun
    use utility_fun
  
  implicit none

  integer, dimension(2) :: numTests, numAsserts, numAssertsTested, numFailures

    write(*,*)
  write(*,*) "combinatory test suite:"
  call test_combinatory &
    ( numTests(1), numAsserts(1), numAssertsTested(1), numFailures(1) )
  write(*,1) numAssertsTested(1), numAsserts(1), &
    numTests(1)-numFailures(1), numTests(1)
  1 format('Passed ',i0,' of ',i0,' possible asserts comprising ',i0,' of ',i0,' tests.')
    write(*,*)
  write(*,*) "utility test suite:"
  call test_utility &
    ( numTests(2), numAsserts(2), numAssertsTested(2), numFailures(2) )
  write(*,1) numAssertsTested(2), numAsserts(2), &
    numTests(2)-numFailures(2), numTests(2)
  2 format('Passed ',i0,' of ',i0,' possible asserts comprising ',i0,' of ',i0,' tests.')
  
  write(*,*)
  write(*,'(a)') "==========[ SUMMARY ]=========="
      write(*,'(a9)',advance="no") " combinatory:"
  if ( numFailures(1) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,'(a9)',advance="no") " utility:"
  if ( numFailures(2) == 0 ) then
    write(*,*) " passed"
  else
    write(*,*) " failed   <<<<<"
  end if
    write(*,*)

  if ( sum(numFailures) /= 0 ) stop 1

end program TestRunner
