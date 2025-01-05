program montecarlo

	!* **********************************************************************!
	!*              Diffusion Monte Carlo                                   *!
	!*                 Sample Program                                       *!
	!*                                                                      *!
	!*  This program pefforms Monte Carlo Integration on the Schrodinger 	*!
	!*  wave equation for several systems, including the H.O., LJ. P., 	*!
	!*  and for H, H2+, and H2.						*!
	!*                                                                      *!
	!*   (c) 1995, Byron Faber, Ioan Kosztin and Klaus Schulten             *!
	!*             University of Illinois at Urbana Champaign		*!
	!*   (Fortran) 2024, Erick Gabriel Fernandes Farias			*!
	!*                   Universidade Federal de Goiás    			*!           
	!* **********************************************************************!

use dmc_module

   	!* **************************************************************!
   	!* VARIABLE DEFINITIONS						*!
   	!*								*!
   	!* **psips = replica matrix					*!
   	!* *ener = array used to store Eo values (for deviation)	*!
   	!* *wz = `boxes' used to sort replicas				*!
   	!* itime = measures time					*!
   	!* rf & ri = right & left measuring bounds for replica boxes 	*!
   	!* others = misc                                             	*!
   	!* **************************************************************!

integer(8) :: potnum
real*8, allocatable :: psips(:,:), wz(:), ener(:)
integer(8) :: itime
real*8 :: vrefav, deviation, ener0, dev0, tmcount
integer :: count_rate, count1, count2
real :: elapsed_time

	!**********      USER initialization    ***********!

call system("clear")
write(*, *) "*****************************************************"
write(*, *) "*                                                   *"
write(*, *) "*              Diffusion Monte Carlo                *"
write(*, *) "*                   Simulator                       *"
write(*, *) "*                                                   *"
write(*, *) "*****************************************************"
write(*, *)
write(*, *) "A.  one-dimensional potentials"
write(*, *) "  1.  Harmonic Oscillator"
write(*, *) "  2.  Morse Potential"
write(*, *) "------------------------------"
write(*, *) "B.  3-dimensional potentials"
write(*, *) "  3.  Hydrogen atom"
write(*, *) "  4.  H2+ molecule"
write(*, *) "  5.  H2 molecule"
write(*, *) "------------------------------"
write(*, *) "Select a potential to simulate (1-5): "
do
	read(*,*) potnum        
	if (potnum >= 1 .and. potnum < 6) then
		exit
	else
		write (*, *) "Invalid value, try again."
	end if 	
end do

if (potnum < 3) then
	Ndim = 1
	else if (potnum > 2 .and. potnum < 5) then
        	Ndim = 3
	else if (potnum == 5) then
		Ndim = 6
	else
		print *, "Error."
end if

write(*, *) "Number in parenthesis are suggestions"
write(*, *) "Ndim = ", Ndim
write(*, *) "Enter Number of Replicas (500):"
read(*, *) Npsips
write(*, *) "Replicas = ", Npsips
write(*, *) "Enter the max number of replicas (2000):"
read(*, *) Nmax
write(*, *) "Nmax = ", Nmax
write(*, *) "Enter random number seed:"
read(*, *) seed
if (seed > 0) seed = -seed
write(*, *) "Seed = ", seed
write(*, *) "Enter the amount of time to run the simulation (1000):"
read(*, *) Nchec
if (potnum .ne. 3) then
	write(*, *) "Enter the left limit for sampling data (-20):"
	read(*, *) ri
else
        ri = 0.0
end if
write(*, *) "left = ", ri
write(*, *) "Enter the right limit for sampling data (20):"
read(*, *) mrf
write(*, *) "right = ", mrf
write(*, *) "Enter the number of boxes to sort into (200):"
read(*, *) Nwf
Ncol = Nchec
write(*, *) "boxes = ", Nwf
write(*, *) "Enter a time step (.1):"
read(*, *) dt
	
if (potnum < 6 .and. potnum > 3) then
	write(*, *) "Enter deviation from R (H - H distance) in units of bohr radius:"
	read(*, *) deviation
	if (potnum == 4) then
		rhh = rhh2 + deviation
	else if (potnum == 5) then
		rhh = crhh + deviation
	end if
	write(*, *) "Using ", rhh, " as hydrogen-hydrogen distance."
end if

call SYSTEM_CLOCK(count_rate=count_rate)
if (count_rate == 0) then
	print *, "Error: SYSTEM_CLOCK does not support accurate counting."
        stop
end if
call SYSTEM_CLOCK(count1)

	! Initialize all arrays with user's entered data   
    
ener = dvector(1_8, Ncol)
psips = dmatrix(1_8 , Nmax, 1_8 , 1 + Ndim)
wz = dvector(1_8, Nwf)

	! Open all files for output

open(unit=1, file="Eo.dat", status="unknown", action="write")
open(unit=2, file="waveshape.dat", status="unknown", action="write")
open(unit=3, file="wavedump.dat", status="unknown", action="write")
write(1, *) "# Format:"
write(1, *) "# column 1:  time dimensionless units"
write(1, *) "# column 2:  energy (unitless for H.O., M.P., Hydrogen, eV otherwise)"
write(2, *) "# Format:"
write(2, *) "# column 1:  position (unitless for H.O, M.P., units of Bohr radius otherwise)"
write(2, *) "# column 2:  R (for hydrogen), in other cases the wave function normalized"
write(2, *) "# column 3:  r * R (for hydrogen), does NOT occur otherwise"
write(2, *) "# Note:  For hydrogen we show the radial wavefunction, all other cases show a distribution"
write(2, *) "# across a single axis"
write(3, *) "# this is a dump of the replica matrix at the end of the simulation."
write(3, *) "# x, y, z for all electrons"
write(1, *) "# psips = ", Npsips, " seed = ", seed, " dt = ", dt, " time = ", Nchec, " left = ", ri
write(1, *) "# right = ", mrf, " boxes = ", Nwf
write(2, *) "# psips = ", Npsips, " seed = ", seed, " dt = ", dt, " time = ", Nchec, " left = ", ri
write(2, *) "# right = ", mrf, " boxes = ", Nwf
write(3, *) "# psips = ", Npsips, " seed = ", seed, " dt = ", dt, " time = ", Nchec, " left = ", ri
write(3, *) "# right = ", mrf, " boxes = ", Nwf
vref = 0.0_8

write(*, *) "Psips initialized"
call initpsips(psips, Nmax, 1+Ndim, potnum, Npsips, Nmax, rhh)

	!**********************************************************************!
	!** Step # 1:  Converge on Eo                                        **!
	!**                                                                  **!
	!** This is the first step in our Monte Carlo method.  Here we       **!
	!** move our replicas, which were created in psips[m][n], by calling **!
	!** a function called walk().  We then kill them or create them      **!
	!** using a function called branch().  At the end of this process we **!
	!** will have converged to our ground state energy Eo, assuming that **!
	!** we have waited a long enough time.                               **!
	!**                                                                  **!
	!**********************************************************************!

write(*, *) "Step 1: Converge to Eo"
itime = 0
vrefav = 0.0
do i = 1, Nchec
	call walk(psips, Nmax, 1+Ndim, potnum, dt, seed)
	call branch(psips,Nmax, 1+Ndim , potnum, Npsips, seed, dt, Ndim, rhh, vref)
	itime = itime + 1
	vrefav = vrefav + vref
	if (potnum > 3) then
		write(1, *) itime * dt, (vrefav / i) * harttoev
	else
		write(1, *) itime * dt, vrefav / i
    	end if
end do
vrefav = vrefav / Nchec

	!*******************************************************************!
	!**  Step 2:  Find wave function                                  **!
	!**                                                               **!
	!** Now that we have converged to our ground state energy we      **!
	!** can count the replicas and where they are.  This can be done  **!
	!** easily by summing the positions of the replicas over a        **!
	!** function of time.  This is more efficient then moving a large **!
	!** population of replicas because at the end of the summation    **!
	!** the wave shape will look the same.  Its a slight trick, but   **!
	!** save a lot of processor time.                                 **!
	!**                                                               **!
	!** Note:  This time we do the same thing as above, walk(), then  **!
	!**        branch(), except this time we also count().  Each of   **!
	!**        these functions is explained in the function itself    **!
	!*******************************************************************!
    
write(*, *) "Step 2: Find wave function"

stwf = (mrf - ri) / (Nwf - 1)

	! Averaging and building wave function
wz = 0.0

do id = 1, Ncol
	call walk(psips, Nmax, 1+Ndim, potnum, dt, seed)
	call branch(psips,Nmax, 1+Ndim , potnum, Npsips, seed, dt, Ndim, rhh, vref)
        call count(wz, psips, Nmax, 1 + Ndim, potnum, ri, Nwf, stwf)
        ener(id) = vref
        itime = itime + 1
end do

	!* Calculate Eo & its deviation *!
  
ener0 = sum(ener) / Ncol
dev0 = sqrt(sum((ener - ener0)**2) / Ncol)

if (potnum > 2) then
	print *, "E0 = ", ener0 * harttoev, "eV    deviation = ", dev0 * harttoev
end if

if (potnum < 4) then
	print *, "E0 = ", ener0, "    deviation = ", dev0
end if

	! Re-normalize
tmcount = 0.0d0

  if (potnum < 3 .or. potnum > 3) then
    do ij = 1, Nwf
      tmcount = tmcount + wz(ij) * wz(ij)
    end do
    
    do ij = 1, Nwf
      write(2, *) (ij - 1) * stwf + ri, wz(ij) / sqrt(tmcount * stwf)
    end do
  end if

  if (potnum == 3) then
    do ij = 1, Nwf
      wz(ij) = wz(ij) / sqrt(4.0d0 * PI * (ij * stwf)**2)
    end do
    
    do ij = 1, Nwf
      tmcount = tmcount + (wz(ij)**2)
    end do
    do ij = 1, Nwf
      write(2, *) (ij - 1) * stwf + ri, wz(ij) / (sqrt(tmcount * stwf) * ij * stwf), wz(ij) / sqrt(tmcount * stwf)
    end do
  end if

  if (potnum > 2) then
    do i = 1, Nmax
      if (psips(i, 1) == 1.0d0) then
        write(3, *) psips(i, 2), psips(i, 3), psips(i, 4)
      end if
    end do
  end if

  if (potnum > 4) then
    do i = 1, Nmax
      if (psips(i, 1) == 1.0d0) then
        write(3, *) psips(i, 5), psips(i, 6), psips(i, 7)
      end if
    end do
  end if

deallocate(wz)
deallocate(ener)
deallocate(psips)
	close(unit=1)
	close(unit=2)
	close(unit=3)

call SYSTEM_CLOCK(count2)
end_time = real(count1) / real(count_rate)
elapsed_time = real(count2 - count1) / real(count_rate)
print *, "Total time of aplication:", elapsed_time, "s"

end program montecarlo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SUBROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!***********************************************************************!
	!**  Intialize the replica matrix:                                    **!
	!**                                                                   **!
	!**  This function sets up the replica matrix that we created in our  **!
	!**  main function.  Simply, we have a Nmax by (1 + Ndim) matrix,     **!
	!**  where Nmax is the maximum number of replicas and (1 + Ndim) is   **!
	!**  the number of dimensions that our replicas `live' in.            **!
	!**                                                                   **!
	!**  Here we have several cases:                                      **!
	!**  harmonic oscillator:  replicas start at 0 on x-axis              **!
	!**  Morse Potential:      replicas start at 0                        **!
	!**  H atom:               replicas start at 1 bohr radius            **!
	!**  H2+ atom:             replicas start at 1 bohr radius            **!
	!**  H2  atom:             replicas start at 2x rhh distance, which   **!
	!**                        is the distance separating the two H       **!
	!**                        molecules in H2.                           **!
	!**                                                                   **!
	!**  Note:  the first two problems are worked in dimensionless units, **!
	!**         while the last ones are worked in units of the bohr radius**!
	!**         The answers for Eo given in the first case are dimension- **!
	!**         less, while the ones in the last case are converted to eV **!
	!**                                                                   **!
	!***********************************************************************!

subroutine initpsips(psips, m, n, potnum, Npsips, Nmax, rhh)
implicit none

integer(8), intent(in) :: m, n, potnum, Npsips, Nmax
real(8), allocatable, intent(inout) :: psips(:,:)
real(8), intent(in) :: rhh
integer(8) :: i

	!* -------------> Initial distribution of Psips */
if (potnum == 5) then
	do i = 1, Npsips
		psips(i, 1) = 1.0
		psips(i, 2) = 0.0
		psips(i, 3) = 0.0
		psips(i, 4) = rhh
		psips(i, 5) = 0.0
		psips(i, 6) = 0.0
		psips(i, 7) = -rhh
	end do

	do i = Npsips + 1, Nmax
		psips(i, 1) = 0.0
        end do
end if

if (potnum < 3) then
	do i = 1, Npsips
		psips(i, 1) = 1.0
		psips(i, 2) = 0.0
        end do

        do i = Npsips + 1, Nmax
		psips(i, 1) = 0.0
		psips(i, 2) = 0.0
        end do
end if

if (potnum < 5 .and. potnum > 2) then
	do i = 1, Npsips
        	psips(i, 1) = 1.0
        	psips(i, 2) = 0.0
        	psips(i, 3) = 0.0
        	psips(i, 4) = 1.0
	end do

	do i = Npsips + 1, Nmax
        	psips(i, 1) = 0.0
	end do
end if

end subroutine initpsips

	!*********************************************************************!
	!** Counting algorithm                                              **!
	!** This is the algorithm used to count replicas & plot our wave-   **!
	!** function.  In order to count them we first take the position of **!
	!** each replica, subtract it from our left boundry, and then       **!
	!** divide by the distance per box.  We can then use this number as **!
	!** an index for an array that we are sorting our replicas into.    **!
	!**                                                                 **!
	!** When reporting our positions back to the user we reverse the    **!
	!** process.                                                        **!
	!*********************************************************************!

subroutine count(wz, psips, m, n, potnum, ri, Nwf, stwf)

use nrutil_module

implicit none
    
integer(8), intent(in) :: m, n, potnum, Nwf
real(8), intent(in) :: ri, stwf
real(8), allocatable, intent(inout) :: psips(:,:), wz(:)
integer :: ij, iz1
real(8) :: r

	! Count the replicas and make a histogram
if (potnum < 3 .or. potnum > 3) then
	do ij = 1, m
		if (psips(ij, 1) == 1.0d0) then
			iz1 = int((psips(ij, 2) - ri) / stwf) + 1
			if (iz1 < Nwf .and. iz1 > 1) then
         			wz(iz1) = wz(iz1) + 1.0d0
			end if
       		end if
	end do
end if

	! For the Hydrogen atom we sample the radius!!!
if (potnum == 3) then
	do ij = 1, m
		if (psips(ij, 1) == 1.0d0) then
                	r = sqrt(DSQR(psips(ij, 2)) + DSQR(psips(ij, 3)) + DSQR(psips(ij, 4)))
                	iz1 = (r / stwf) + 1_8
                	if (iz1 >= 1 .and. iz1 <= Nwf) then
                    		wz(iz1) = wz(iz1) + 1.0d0
                	end if
            	end if
        end do
end if

end subroutine count

	!*********************************************************!
	!** Walk:  guassian movement                            **!
	!**                                                     **!
	!** This function moves the replica distribution by a   **!
	!** guassian distributed randomn number.                **!
	!**                                                     **!
	!** i.e:  Each variable gets moved by a value           **!
	!**       sqrt(dt) * guassian number from its original  **!
	!**       location.                                     **!
	!**                                                     **!
	!**  Note that we don't have to worry about units       **!
	!**  because we are working in unitless coordinates for **!
	!**  all of our cases.  For 1-d cases we report unitless**!
	!**  answers.  For 3-d we convert our values to eV at   **!
	!**  the end.                                           **!
	!*********************************************************!

subroutine walk(psips, m, n, potnum, dt, seed)

use gasdev_module, only: gasdev

implicit none

integer(8), intent(in) :: m, n, potnum
integer(8), intent(inout) :: seed
real(8), intent(in) :: dt
real(8), allocatable, intent(inout) :: psips(:,:)
integer(8) :: k, ii
real(8) :: sigma

sigma = sqrt(dt)

do k = 1, m
	if (psips(k, 1) == 1.0d0) then
        	do ii = 1, n - 1 
			if (ii + 1 <= n) then
				psips(k, ii + 1) = psips(k, ii + 1) + sigma * gasdev(seed)
			else
				print *, "Error: ii + 1 goes beyond the limits of psips matrix!"
				stop
			end if
		end do
	end if
end do

end subroutine walk

	!****************************************************************!
	!** Branching                                                  **!
	!**                                                            **!
	!** This function takes our replica matrix and does the        **!
	!** branching process on it.  This is done in the following    **!
	!** order:                                                     **!
	!** 1.)  Calculate the Vlocal and Vaverage of the replicas by  **!
	!**      using the right potential.                            **!
	!** 2.)  Use the Vlocal and compare it to Vaverage to find out **!
	!**      if the replica will be born or will die.  If the      **!
	!**      replica dies, we simply delete it from the matrix.    **!
	!**      However, if a new one is born, then we must search    **!
	!**      the matrix for an empty spot, and then enter the      **!
	!**      new replica in.  (This replica will have the coord.   **!
	!**      of its parent particle).                              **!
	!****************************************************************!

subroutine branch(psips, m, n, potnum, Npsips, seed, dt, Ndim, rhh, vref)

	!* Now, find all the replicas, count them using nact & find the potential *!

use ran3_module, only: ran3
use potentials_module, only: vho, vmp, vh, vh2p, vh2
use nrutil_module
  
implicit none

integer(8), intent(in) :: m, n, potnum, Npsips, Ndim
integer(8), intent(inout) :: seed
real(8), intent(in) :: dt, rhh
real(8), intent(inout) :: vref
real(8), allocatable, intent(inout) :: psips(:,:)
real(8), allocatable :: vtmp(:), cor(:)
integer(8) :: kdim, inew, nnewp, ja, ib, jj, npb, ii, nact
real(8) :: vloc, pb, fpb, vav
logical :: tobir

vav = 0.0
nact = 0
vref = 0.0

vtmp = dvector(1_8, m)
cor = dvector(1_8, n - 1_8)
  
do ja = 1, m
	if (psips(ja, 1) == 1.0) then
		nact = nact + 1
		cor(:) = psips(ja, 2:n)
		select case (potnum)
			case (1)
        			vtmp(ja) = vho(cor, n - 1_8)
      			case (2)
        			vtmp(ja) = vmp(cor, n - 1_8)
      			case (3)
        			vtmp(ja) = vh(cor, n - 1_8, Ndim)
      			case (4)
        			vtmp(ja) = vh2p(cor, n - 1_8, Ndim, rhh)
      			case (5)
        			vtmp(ja) = vh2(cor, n - 1_8, rhh)
      			case default
        			stop
      		end select
      		vav = vav + vtmp(ja)
	end if
end do

if (nact > 0) then
	vav = vav / nact
end if

vref = vav - (nact - Npsips) / (Npsips * dt)

do jj = 1, m
	if (psips(jj, 1) == 1.0d0) then
		vloc = vtmp(jj)

      		if (vloc < vref) then
	!* -----------------> birth
      			pb = exp(-(vloc - vref) * dt) - 1.0d0
        		if (pb >= 2.0d0) pb = 2.0d0
        			npb = floor(pb)
        			fpb = pb - npb

        			if (fpb >= ran3(seed)) then
          				nnewp = npb + 1
        			else
          				nnewp = npb
        			end if

        ! Loop to create new replicas
			do inew = 1, nnewp
          			ib = 1
          			tobir = .true.
          			do while (tobir)
            				if (psips(ib, 1) == 1.0d0) then
              				ib = ib + 1
              				if (ib > m) then
                				print *, "Maximum number of replicas exceeded"
                				stop
              				end if
            				else
              					psips(ib, 1) = 1.0d0
              					psips(ib, 2:n) = psips(jj, 2:n)
              					vtmp(ib) = vref
              					tobir = .false.
            				end if
          			end do
			end do
     			else
	!* -----> dead
				pb = 1.0d0 - exp(-(vloc - vref) * dt)

        				if (pb > ran3(seed)) then
          					psips(jj, 1) = 0.0d0
        				end if
      		end if
	end if
end do

call free_dvector(vtmp, 1_8, m)
call free_dvector(cor, 1_8, n - 1_8)

end subroutine branch

