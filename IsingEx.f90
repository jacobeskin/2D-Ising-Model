program IsingEx
  !compilation: make -f makeising
  !running: ./IsingEx.exe

  !Simulation of 2D Ising model with Metropolis Monte Carlo method. Periodic
  !boundary conditions, a site is affected by 4 nearest neighbours. Interaction
  !energy J is to be changed from the code itself, and J=1 is ferromagnet and
  !J=-1 is antiferromagnet. Temperature T and magnetic field strength H are
  !also to be changed if need be straight form the code.Lattice is 100x100
  !sites. Outputs .xyz files that are readable with e.g. Ovito for making
  !visualizations on how the magnetization changes in different regions of the
  !lattice. Also outputs files for plotting changes in magnetization as a
  !function of temperature and magnetic field strength. 

  
  !Used modules
  use mtdefs
  use mtmod
  use IsingMod
  implicit none

  double precision :: J, H, T, start, end
  double precision :: M_ave1, M_ave2, M_ave3, M_ave4, M_ave5, M_ave6
  integer :: i, u, movie, N_MC, lattice(100,100)
  character(len=80) :: filename

  

  ! Initialize RNG
  call sgrnd(getseed(info=1))
  
  !------------------------------------
  ! Decide what J, H and T you want
  !------------------------------------
  J = 1.0d+0
  call cpu_time(start)
  ! Calculate <M> and stuff
  H = 0.2d+0
  T = 0.1d+0
  N_MC = 1000000
  movie = 0

  open(unit=1, file='M_aveVSTH_3', status='unknown')
  open(unit=2, file='M_aveVSTH_4', status='unknown')
  do i = 1,20
     call gen_lattice(lattice)
     call MMC(lattice, 1.0d+0, 0.0d+0, T, N_MC, M_ave1, 1, 'case1.xyz')
     call MMC(lattice, 1.0d+0, H, 0.1d+0, N_MC, M_ave2, movie, 'case2.xyz')
     call MMC(lattice, 1.0d+0, H, 5.0d+0, N_MC, M_ave3, movie, 'case3.xyz')
     call MMC(lattice, -1.0d+0, 0.0d+0, T, N_MC, M_ave4, 1, 'case4.xyz')
     call MMC(lattice, -1.0d+0, H, 0.1d+0, N_MC, M_ave5, movie, 'case5.xyz')
     call MMC(lattice, -1.0d+0, H, 5.0d+0, N_MC, M_ave6, movie, 'case6.xyz')
     M_ave1 = abs(M_ave1)
     M_ave2 = abs(M_ave2)
     M_ave3 = abs(M_ave3)
     M_ave4 = abs(M_ave4)
     M_ave5 = abs(M_ave5)
     M_ave6 = abs(M_ave6)
     write(1, '(f10.5, 4x, f10.5, 4x, f10.5, 4x, f10.5, 4x, f10.5)'), T,&
          &M_ave1, H, M_ave2, M_ave3 
     write(2, '(f10.5, 4x, f10.5, 4x, f10.5, 4x, f10.5, 4x, f10.5)'), T,&
          &M_ave4, H, M_ave5, M_ave6 
     T = T+0.5d+0
     H = H+0.2d+0
     movie = 1
  end do
  close(1)
  close(2)
  
  call cpu_time(end)

  print*, (end-start)/60 


end program IsingEx
     
       
       
             
          

       

       

       

    
    
    
