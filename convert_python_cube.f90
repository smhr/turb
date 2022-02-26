      module conversion
      contains
        subroutine convert(fin, fout)
        use iso_c_binding
        implicit none
        
        integer, parameter            :: DPR = selected_real_kind(p=15)
        character(len=*), intent(in)  :: fin
        character(len=*), intent(in)  :: fout
        real(DPR),    allocatable, dimension(:, :, :) :: vgrid
        real(kind=4), allocatable, dimension(:, :, :) :: vgrid_s
        integer, parameter :: n = 32
        integer            :: i, j, k
        
        i = 0; j = 0; k = 0
        allocate(vgrid(n, n, n))
        open(10, file=fin, status='old', &
            access='stream', form='unformatted')
        print*,fin
        read(10) vgrid
        close(10)
        
        open(20, file=fout, form='unformatted')
        vgrid_s = vgrid ! convert to single precision
        write(20) (((vgrid_s(i,j,k),i=1,n),j=1,n),k=1,n)
        close(20)
        
        end subroutine convert
      end module conversion

!================================
      program convert_python_cube
      use conversion
      implicit none

      call convert('vx.bin', 'cube_v1.dat')
      call convert('vy.bin', 'cube_v2.dat')
      call convert('vz.bin', 'cube_v3.dat')


      end program convert_python_cube
