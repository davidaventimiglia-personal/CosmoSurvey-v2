module utils
  implicit none
contains
  ! Output data appropriate for a GNUPlot 3-D histogram
  subroutine hist3d(dist, x1, x2, y1, y2, file, floor)
    real, dimension(:, :), intent(in) :: dist
    real, intent(in) :: x1, x2, y1, y2
    integer, intent(in) :: file
    real, optional, intent(in) :: floor
    integer :: i, j, M, N
    real :: f
    real, dimension(size(dist, 2)) :: floorline
    f = 1e-10
    if (present(floor)) f = floor
    M = size(dist, 1)
    N = size(dist, 2)
    floorline = [(f, i = 1, N)]
    do i = 1, M
       call isoline(floorline, x1, x2, y1, y2, file, i-1, f)
       call isoline(dist(i, :), x1, x2, y1, y2, file, i-1, f)
       call isoline(dist(i, :), x1, x2, y1, y2, file, i, f)
       call isoline(floorline, x1, x2, y1, y2, file, i, f)
    end do
  end subroutine hist3d

  ! Output data for a GNUPlot isoline
  subroutine isoline(iso, x1, x2, y1, y2, file, i, floor)
    real, dimension(:), intent(in) :: iso
    real, intent(in) :: x1, x2, y1, y2
    integer, intent(in) :: file
    integer, intent(in) :: i
    real, intent(in) :: floor
    real :: dx, dy
    integer :: j, N
    N = size(iso)
    dx = (x2 - x1)/real(N)
    dy = (y2 - y1)/real(N)
    do j = 1, N
       write(file, *) x1 + real(i)*dx, y1 + real(j-1)*dy, floor
       write(file, *) x1 + real(i)*dx, y1 + real(j-1)*dy, iso(j)
       write(file, *) x1 + real(i)*dx, y1 + real(j)*dy, iso(j)
       write(file, *) x1 + real(i)*dx, y1 + real(j)*dy, floor
    end do
    write(file, *) ''
  end subroutine isoline

  ! Regrid data into "bins" represented by truncated pyramids of
  ! negligibly-sloping sides, appropriate for use by PLplot.
  subroutine hist_3d(x, y, z)
    real, intent(inout), dimension(:), allocatable :: x, y
    real, intent(inout), dimension(:, :), allocatable :: z
    real, dimension(:), allocatable :: temp_x, temp_y
    real, dimension(:, :), allocatable :: temp_z
    integer :: Nx, Ny, i, j
    Nx = size(z, 1)
    Ny = size(z, 2)
    allocate(temp_x(4*Nx))
    allocate(temp_y(4*Ny))
    allocate(temp_z(4*Nx, 4*Ny))
    temp_z = minval(z)
    do i = 1, Nx
       do j = 1, Ny
          temp_z(4*(i-1)+2, 4*(j-1)+2) = z(i,j)
          temp_z(4*(i-1)+2, 4*(j-1)+3) = z(i,j)
          temp_z(4*(i-1)+3, 4*(j-1)+2) = z(i,j)
          temp_z(4*(i-1)+3, 4*(j-1)+3) = z(i,j)
       end do
    end do
    do i = 1, Nx
       temp_x(4*(i-1)+1) = x(i)-((x(Nx)-x(1))/real(Nx-1))/2.*0.999
       temp_x(4*(i-1)+2) = x(i)-((x(Nx)-x(1))/real(Nx-1))/2.*0.998
       temp_x(4*(i-1)+3) = x(i)+((x(Nx)-x(1))/real(Nx-1))/2.*0.998
       temp_x(4*(i-1)+4) = x(i)+((x(Nx)-x(1))/real(Nx-1))/2.*0.999
    end do
    do j = 1, Ny
       temp_y(4*(j-1)+1) = y(j)-((y(Ny)-y(1))/real(Ny-1))/2.*0.999
       temp_y(4*(j-1)+2) = y(j)-((y(Ny)-y(1))/real(Ny-1))/2.*0.998
       temp_y(4*(j-1)+3) = y(j)+((y(Ny)-y(1))/real(Ny-1))/2.*0.998
       temp_y(4*(j-1)+4) = y(j)+((y(Ny)-y(1))/real(Ny-1))/2.*0.999
    end do
    deallocate(x)
    deallocate(y)
    deallocate(z)
    allocate(x(4*Nx))
    allocate(y(4*Ny))
    allocate(z(4*Nx, 4*Ny))
    x = temp_x
    y = temp_y
    z = temp_z
  end subroutine hist_3d

  character(len=80) function f_nvpair (name)
    character(len=*) :: name
    character(len=80) :: arg
    integer :: i
    f_nvpair = 'NOT FOUND'
    do i = 1, iargc()
       call getarg(i, arg)
       if ((arg.eq.trim(name)).and.(iargc().gt.i)) call getarg(i+1, f_nvpair)
    end do
  end function f_nvpair

  function f_split (string, char)
    character(len=*), intent(in) :: string
    character(len=1), intent(in) :: char
    character(len=80), dimension(2) :: f_split
    character(len=80) :: a, b
    integer :: idx
    a = ''
    b = ''
    idx = index(string, char)
    if (idx.gt.0) then
       if (idx.gt.1) then
          a = string(1:idx-1)
       end if
       if (idx.lt.80) then
          b = string(idx+1:80)
       end if
    end if
    f_split(1) = trim(a)
    f_split(2) = trim(b)
  end function f_split

  function f_records (filename)
    character(len=*), intent(in) :: filename
    character(len=80), dimension(:), pointer :: f_records
    character(len=80) :: record
    integer :: N, i
    open (unit=10, file=filename, action='read')
    N = f_record_count(10)
    allocate(f_records(N))
    do i = 1, N
       read (unit=10, fmt='(A)') record
       f_records(i) = record
    end do
    close (unit=10)
  end function f_records

  integer function f_record_count (fu)
    integer, intent(in) :: fu
    character(len=80) :: rec
    integer :: stat, i
    stat = 0
    f_record_count = 0
    rewind (unit = fu)
    do while (stat.ge.0)
       read (unit=fu, fmt='(A)', iostat=stat) rec
       if (stat.ge.0) f_record_count = f_record_count + 1
    end do
    rewind (unit = fu)
  end function f_record_count

  ! Output matrix m to unit in JSon format
  subroutine json_matrix (m, unit)
    real, dimension(:) :: m
    integer :: i, j, unit
    write(unit, '(a)', advance='no') '['
    do i = 1, size(m, 1) - 1
       write(unit, '(e10.3)', advance='no') m(i)
       write(unit, '(a)', advance='no') ','
    end do
    write(unit, '(e10.3)', advance='no') m(size(m, 1))
    write(unit, '(a)') ']'
  end subroutine json_matrix
end module utils
