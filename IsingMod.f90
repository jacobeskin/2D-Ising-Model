module IsingMod
use mtdefs
use mtmod

  implicit none

  ! Parameters
  double precision, parameter :: kb = 0.00008613324

  contains
  
  !------------------------------------- 
  ! Subroutine for creating the lattice-
  !-------------------------------------
    
    subroutine gen_lattice(lattice)
      implicit none
      
      integer, intent(out) :: lattice(100,100)
      integer :: u, i, j

      do i=1,100
         do j=1,100
            u = igrnd(0,1)
            if (u==0) then 
               u=-1
            end if
            lattice(j,i) = u
         end do
      end do
    end subroutine gen_lattice
  !------------------------------------------

  !------------------------------------------
  ! Subroutine for calculating magnetization-
  !------------------------------------------

    subroutine magnetization(lattice, M)
      implicit none
      
      integer, intent(in) :: lattice(100,100)
      double precision, intent(out) :: M
      double precision :: s_tot
      integer :: i, j, si, N_spins = 10000

      ! Calculate magnetization per spin, M = sum(s_i)/N_spins
      s_tot = 0.0
      do i = 1,100
         do j = 1,100
            si = lattice(j,i)
            s_tot = s_tot+si
         end do
      end do
      M = s_tot/N_spins
      
    end subroutine magnetization
  !-------------------------------------------

  !------------------------------
  ! Subroutine for Metropolis MC-
  !------------------------------
  
    subroutine MMC(lattice, J, H, T, N_MC, M_ave, movie, filename)
      implicit none
      
      integer, intent(in) :: lattice(100,100), N_MC
      integer, intent(in) :: movie ! Parameter to decide if we want visualize
      double precision, intent(in) :: H, T, J
      
      double precision, intent(out) :: M_ave
      
      double precision :: M, M_sum, dE, u
      integer :: i, l, k, x, y, s, sn, x1, x2, y1, y2, lat(100,100)
      character(len=9), intent(in) :: filename

      lat = lattice
      
      ! If movie==0, then record this simulation
      if (movie==0) then
         open(10, file=filename, status='unknown')
      end if
      
      ! Metropolis routine
      M_sum = 0.0
      do i = 1,N_MC
         
         ! Choose lattice site randomly
         x = igrnd(1,100)
         y = igrnd(1,100)
         
         s = lat(x,y)
         sn = -1*s
         
         ! Assigning nearest neighbors, using periodic boundary conditions if 
         ! necessary
         
         ! Neighbor above:
         if (x==1) then ! If we are on the firs row of the lattice
            y1 = lat(100,y)
         else           
            y1 = lat(x-1,y)
         end if
         
         ! Neighbor below:
         if (x==100) then ! If we are on the last row of the lattice
            y2 = lat(1,y)
         else
            y2 = lat(x+1,y)
         end if
         
         ! Neighbor to the left:
         if (y==1) then ! If we are on the first column of the lattice
            x1 = lat(x, 100)
         else
            x1 = lat(x,y-1)
         end if
         
         ! neihbor to the right:
         if (y==100) then ! If we are on the last column of the lattice
            x2 = lat(x,1)
         else
            x2 = lat(x,y+1)
         end if

         ! Calculate the change in energy
         dE = -2.0*J*(x1+x2+y1+y2)*sn-H*sn
         
         ! Deciding if we keep the flip
         if (dE<=0) then
            lat(x,y) = sn
         else
            u = grnd()
            if (u<exp(-1.0*dE/T)) then
               lat(x,y) = sn
            end if
         end if
         
         ! Calculate Magnetization
         call magnetization(lat, M)
         
         ! Update M_sum
         M_sum = M_sum+M
         
         ! Print out xyz for visualization
         if (movie==0) then
            if (mod(i,10000)==1) then
               write(10,*) 10000
               write(10,*) 'Step', i
               do l=1,100
                  do k=1,100
                     if (lat(k,l)==1) then
                        write(10,*) "U",l,k
                     else
                        write(10,*) "D",l,k
                     end if
                  end do
               end do
               call flush(10)
            end if
         end if
         
      end do
      if (movie==0) then
         close(10)
      end if
      
      ! Calculate average of magnetizations per spin
      M_ave = M_sum/N_MC
      
    end subroutine MMC
    
  end module IsingMod
