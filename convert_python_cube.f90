! This program converts make_turb output files (vx.bin, vy.bin, vz.bin)
!to a format that phantom can read.
! Please note that here n must be equal to half of res in make_turb.py script.

      module conversion
      contains
        subroutine convert(fin, fout)
        use iso_c_binding
        use iso_fortran_env, only : file_storage_size
        implicit none
        
        integer, parameter            :: DPR = selected_real_kind(p=15)
        character(len=*), intent(in)  :: fin
        character(len=*), intent(in)  :: fout
        real(DPR),    allocatable, dimension(:, :, :) :: vgrid
        real(kind=4), allocatable, dimension(:, :, :) :: vgrid_s
        integer            :: res
        integer            :: i, j, k
        integer            :: file_size
        
        i = 0; j = 0; k = 0
        inquire(file=fin, size=file_size)
        res = nint((file_size * file_storage_size /8)**(1.0/3.0))/2
        allocate(vgrid(res, res, res))
        open(10, file=fin, status='old', &
            access='stream', form='unformatted')
        
        write(*,*)'size of file '//fin//' is ',file_size * FILE_STORAGE_SIZE /8,' bytes'
        
        print"(1x,a,2i5)", "Resolution =", res

        read(10) vgrid
        close(10)
        
        open(20, file=fout, form='unformatted')
        vgrid_s = vgrid ! convert to single precision
        write(20) (((vgrid_s(i,j,k),i=1,res),j=1,res),k=1,res)
        close(20)
        
        end subroutine convert
      end module conversion

!================================
      program convert_python_cube
      use conversion
      implicit none
      !integer, parameter :: res = 128 ! Cube resolution. E.g. res = 128 means a 128^3 velocity cube. 
                                      ! Must be equal to res in make_turb.py script.

      call convert('vx.bin', 'cube_v1.dat')
      call convert('vy.bin', 'cube_v2.dat')
      call convert('vz.bin', 'cube_v3.dat')
      print"(1x,a)", "Successfully convert all cubes."

      end program convert_python_cube
