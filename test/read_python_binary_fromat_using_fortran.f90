      program read_python_binary_fromat
      use iso_c_binding
      implicit none

      integer, parameter :: DPR = selected_real_kind(p=15)
      character(len=*), parameter :: filename = 'vx.bin'

      integer, parameter :: n = 16
      integer :: i, j, k
      real(DPR), allocatable, dimension(:, :, :) :: a
      
      i = 0; j = 0; k = 0

      allocate(a(n, n, n))

      open(40, file=filename, status='old', &
           access='stream', form='unformatted')

      read(40) a
      close(40)
      do i = 1, n
        do j = 1, n
            write(*,'(f16.8)') (a(i, j, k), k = 1, n)
        enddo
      enddo

      end

