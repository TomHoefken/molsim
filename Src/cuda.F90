module CUDAModule

      use mol_cuda
      use MolModule
      use precision_m
      implicit none

      logical, device :: laccepted

      real(fp_kind), device, allocatable   :: pmetro(:)
      real(fp_kind), device, allocatable   :: ptranx(:)
      real(fp_kind), device, allocatable   :: ptrany(:)
      real(fp_kind), device, allocatable   :: ptranz(:)
      real(fp_kind), device, allocatable   :: pmove(:)
      real(fp_kind), allocatable   :: pmetro_h(:)
      real(fp_kind), allocatable   :: ptranx_h(:)
      real(fp_kind), allocatable   :: ptrany_h(:)
      real(fp_kind), allocatable   :: ptranz_h(:)
      !real(8),device, allocatable    :: dtran_d(:)
     ! real(8),device, allocatable :: prandom
      !integer(4),parameter    :: isamp = 1   ! sequential 0 or random 1
      integer(4),device,allocatable              :: mcstat_d(:,:)     ! 0 for rejected and 1 for accepteD
      integer(4),device              :: imcaccept_d = 1
      integer(4),device              :: imcreject_d = 2
      integer(4),device              :: imcboxreject_d = 3
      integer(4),device              :: imchsreject_d = 4
      integer(4),device              :: imchepreject_d = 5
      integer(4),device              :: imovetype_d
      integer(4), device, allocatable :: ispartmove_d(:)
      !integer(4)              :: ichargechangemove = 2
      integer(4), allocatable :: arrevent(:,:,:)
      real(fp_kind), constant   :: One_d = 1.0d0
      real(fp_kind), constant   :: Zero_d = 0.0d0
      logical  :: lboxoverlap = .false.! =.true. if box overlap
      !logical  :: lhsoverlap  = .false.! =.true. if hard-core overlap
      logical  :: lhepoverlap = .false.! =.true. if hard-external-potential overlap
      real(8)  :: weight      ! nonenergetic weight
    !  real(8),device, allocatable  :: deltaEn(:)
      !real(8)                 :: utot
      !real(8), device         :: utot_d
      integer(4), device :: blocksize = 256
      integer(4)         :: blocksize_h = 256
      integer(4), device             :: lock = 0
      integer(4), device             :: istat
      integer(4), device             :: goal = 5
      integer(4) :: iloops
      integer(4), device :: iloops_d
      integer(4) :: iblock1, iblock2
      integer(4), device :: iblock1_d
      integer(4), device :: iblock2_d
      integer(4) :: ismem

      real(fp_kind), device, allocatable :: Etwo_d(:)
      real(fp_kind), device, allocatable :: Etwoold_d(:)
      real(fp_kind), device, allocatable :: Etwonew_d(:)
      real(fp_kind), device, allocatable :: Ebond_d(:)
      real(fp_kind), device, allocatable :: Eclink_d(:)

      !
      integer(4), device, allocatable :: ipGPU_d(:)

      !MCPass_cuda
      real(fp_kind),device :: dubond_d
      real(fp_kind), device :: duclink_d
      real(fp_kind), device :: dutot_d
      integer(4) :: isharedmem_mcpass

      real(fp_kind),device :: dutwo_d

      !probability for algorithm
      real(8) :: pliang
      real(8) :: pmcpasscuda

      !flag for default MCPass
      logical :: lmcpass = .true.


      contains

      subroutine GPUDriver(iStage)

         implicit none

         integer(4), intent(in) :: iStage
         character(40), parameter :: txroutine ='GPU'
         integer(4) :: istat

         namelist /nmlCUDA/ pliang, pmcpasscuda

         if (ltrace) call WriteTrace(1, txroutine, iStage)

         select case (iStage)
         case (iReadInput)

! ... set initial value
            pliang = 0.0
            pmcpasscuda = 0.0

! ... read input data (nmlParticle)

            rewind(uin)
            read(uin,nmlCUDA)

! ... sum up the probabilities
            pmcpasscuda = pmcpasscuda + pliang


         case (IWriteInput)

! ... allocate memory for arrays on GPU
            call AllocateDeviceParams
            istat = cudaDeviceSynchronize()
            call PrepareMC_cudaAll
            istat = cudaDeviceSynchronize()
            call TransferConstantParams

         case (IBeforeSimulation)

            utot_d = u%tot

         end select


      end subroutine GPUDriver

      subroutine MCPassCUDA

         implicit none
         real(8) :: rrandom, Random
         integer(4) :: ierr, ierra
#if defined (_TESTGPU_)
         real(8) :: Random3

         rrandom = Random3(iseed_test2)
#else
         rrandom = Random(iseed)
#endif

         if (rrandom <= pliang) then
            !call ConvertCoordinatesAuxToDevice<<<iblock1,256>>>
            call MCPassAllGPU
            ierra = cudaDeviceSynchronize()
            call MCCudaUpdate
            lmcpass = .false.
         else if (rrandom < pmcpasscuda) then
            call CalcNewPositions<<<iblock1,256>>>
            call MCPass_cuda<<<iblock2,256,isharedmem_mcpass>>>
            ro = ro_d
            ro_aux = ro
            call ConvertCoordinates<<<iblock1,256>>>
            ierra = cudaDeviceSynchronize()
            call MCCudaUpdate
            lmcpass = .false.
         else
            lmcpass = .true.
         end if


      end subroutine MCPassCUDA


      subroutine MCCudaUpdate

         implicit none
         integer(4) :: ip, ierra
         u%tot = u%tot + dutot_d
         u%twob(0) = u%twob(0) + dutwo_d
         u%bond = u%bond + dubond_d
         u%crosslink = u%crosslink + duclink_d
         ro = ro_aux
         r = ro
         laz = laz_aux_d
         do ip = 1, np
            if (laz(ip)) then
               az(ip) = zat(iatan(ip))
            else
               az(ip) = Zero
            end if
         end do

      end subroutine MCCudaUpdate

         !! subroutine MCPass
         !! main routine in mc step
         !! contains
         !!  subroutines:
         !!              CalcNewPositions
         !!              CalculateUpperPart
         !!              MakeDecision_CalcLowerPart
      subroutine MCPassAllGPU

         use precision_m
         implicit none
         logical, device :: lhsoverlap(np_d)
         integer(4) :: ipart
         integer(4) :: i,j,istat,ierr,ierra

               lhsoverlap = .false.
               Etwo_d = 0.0
               Etwoold_d = 0.0
               !Etwonew_d = 0.0
               Ebond_d = 0.0
               Eclink_d = 0.0
               dutot_d = 0.0
               dutwo_d = 0.0
               dubond_d = 0.0
               duclink_d = 0.0
#if defined (_NORMAL_)
               call GenerateRandoms<<<iblock1,256>>>
#endif
               ierra = cudaDeviceSynchronize()
               call CalcNewPositions<<<iblock1,256>>>
               ierr = cudaGetLastError()
               ierra = cudaDeviceSynchronize()
               if (ierr /= cudaSuccess) write(*,*) "Sync kernel error: ", cudaGetErrorString(ierr)
               if (ierra /= cudaSuccess) write(*,*) "Async kernel err: ", cudaGetErrorString(ierra)
               call CalcWeight<<<iblock1,256>>>
               ierra = cudaDeviceSynchronize()
               call CalculateUpperPart<<<iblock1,256>>>(lhsoverlap)
               ierr = cudaGetLastError()
               ierra = cudaDeviceSynchronize()
               if (ierr /= cudaSuccess) write(*,*) "Sync kernel error: ", cudaGetErrorString(ierr)
               if (ierra /= cudaSuccess) write(*,*) "Async kernel err: ", cudaGetErrorString(ierra)
               do j = 1, iloops
                  ipart = j
                  call MakeDecision_CalcLowerPart<<<iblock2,256>>>(lhsoverlap,ipart,iloops_d)
                  ierr = cudaGetLastError()
                  ierra = cudaDeviceSynchronize()
                  if (ierr /= cudaSuccess) write(*,*) "Sync kernel error: ", cudaGetErrorString(ierr)
                  if (ierra /= cudaSuccess) write(*,*) "Async kernel err: ", cudaGetErrorString(ierra)
                  call CalcLowerPart2<<<iblock2,256>>>(lhsoverlap, ipart, iloops_d)
                  ierr = cudaGetLastError()
                  ierra = cudaDeviceSynchronize()
                  if (ierr /= cudaSuccess) write(*,*) "Sync kernel error: ", cudaGetErrorString(ierr)
                  if (ierra /= cudaSuccess) write(*,*) "Async kernel err: ", cudaGetErrorString(ierra)
               end do
                  ierra = cudaDeviceSynchronize()
               call TransferToAux<<<iblock1,256>>>
                  ierra = cudaDeviceSynchronize()

      end subroutine MCPassAllGPU

      attributes(global) subroutine TransferToAux

         implicit none
         integer(4) :: id

         id = (blockidx%x-1)*blocksize + threadIDx%x
         if (id <= np_d) then
            ro_aux(1,id) = ro_d(1,ipGPU_d(id))
            ro_aux(2,id) = ro_d(2,ipGPU_d(id))
            ro_aux(3,id) = ro_d(3,ipGPU_d(id))
            if (lweakcharge_d) laztm_d(id) = laz_d(id)
            if (lweakcharge_d) laz_aux_d(id) = laz_d(ipGPU_d(id))
         end if

      end subroutine TransferToAux

      subroutine PrepareMC_cudaAll

         use precision_m
         implicit none
         integer(4) :: numblocks
         integer(4) :: sizeofblocks =512



         if (.not.allocated(pmetro)) then
            allocate(pmetro(np))
            pmetro = 0.0
         end if

         if (.not.allocated(ptranx)) then
            allocate(ptranx(np))
            ptranx = 0.0
         end if

         if (.not.allocated(ptrany)) then
            allocate(ptrany(np))
            ptrany = 0.0
         end if

         if (.not.allocated(ptranz)) then
            allocate(ptranz(np))
            ptranz = 0.0
         end if

         if (.not.allocated(pmove)) then
            allocate(pmove(np))
            pmove = 0.0
         end if



         if (.not.allocated(pmetro_h)) then
            allocate(pmetro_h(np))
            pmetro_h = 0.0
         end if

         if (.not.allocated(ptranx_h)) then
            allocate(ptranx_h(np))
            ptranx_h = 0.0
         end if

         if (.not.allocated(ptrany_h)) then
            allocate(ptrany_h(np))
            ptrany_h = 0.0
         end if

         if (.not.allocated(ptranz_h)) then
            allocate(ptranz_h(np))
            ptranz_h = 0.0
         end if

         if (.not.allocated(Etwo_d)) then
            allocate(Etwo_d(np))
            Etwo_d = 0.0
         end if

         if (.not.allocated(Etwoold_d)) then
            allocate(Etwoold_d(np))
            Etwoold_d = 0.0
         end if

         if (.not.allocated(Etwonew_d)) then
            allocate(Etwonew_d(np))
            Etwonew_d = 0.0
         end if

         if (.not.allocated(Ebond_d)) then
            allocate(Ebond_d(np))
            Ebond_d = 0.0
         end if
         if (.not.allocated(Eclink_d)) then
            allocate(Eclink_d(np))
            Eclink_d = 0.0
         end if

         if (.not.allocated(mcstat_d)) then
            allocate(mcstat_d(np,5))
            mcstat_d = 0
         end if

         call CalcWeight_h


         iblock1 = ceiling(real(np) / blocksize_h)
         iloops = ceiling(real(np)/20480)
         iloops_d = iloops
         iblock2 = ceiling(real(iblock1) / iloops)
         iblock1_d = iblock1
         iblock2_d = iblock2

         !ismem = nbuf*fp_kind+fp_kind + fp_kind*npt+npt + fp_kind * nptpt
         !shared memory for MCPass_cuda
           numblocks = floor(Real((nptm*np)/sizeofblocks)) + 1
           isharedmem_mcpass = 2*sizeofblocks*fp_kind + sizeofblocks*4 + threadssum*(nptpt+1)*fp_kind
#if defined (_NORMAL_)
         call GenerateSeeds
#endif

      end subroutine PrepareMC_cudaAll

      !! subroutine GenerateRandoms
      !! Generates the random numbers for decision in Metropolis algorithm and for displacements
      !! contains
      !!  subroutines:
      !!              curandGenerateUnifomr
      !!  internal parameters:
      !!                      gen: test parameter
      !!  module parameters:
      !!                      pmetro(np): probability for acceptance of move(i)
      !!                      ptran(np):  probability for displacement of particles
     attributes(global) subroutine GenerateRandoms

            implicit none
            integer(4) ::  id


            id = (blockidx%x-1)*blocksize + threadIDx%x
               call Random_d(seeds_d(id),pmetro(id),id)
               call Random_d(seeds_d(id),ptranx(id),id)
               call Random_d(seeds_d(id),ptrany(id),id)
               call Random_d(seeds_d(id),ptranz(id),id)
               call Random_d(seeds_d(id),pmove(id),id)



      end subroutine GenerateRandoms


      !! subroutine CalcNewPositions
      !! running on device, calling from device
      !! Calculates the new coordinates for all particles
      !! contains:
      !!  subroutines:
      !!              PBC: calculate coordinates for periodic boundary conditions
      !!              syncthreads: synchronizes all threads
      !!  internal parameters:
      !!              id: global index of thread and index of particle in global list
      !!              blocksize: number of threads in one thread block
      !!  global parameters:
      !!              blockidx%x: Index of thread block on the grid
      !!              threadx%x: internal index of thread in its thread block
      !!              ro(1:3,np): coordinates in x,y,z-direction of the particles
      !!              rotm(1:3,np): new coordinates of the particles
      !!              dtran_d(npt): maximum displacement for particle type ipt
      !!              iptpn(np): particle type of the particles
      attributes(global) subroutine CalcNewPositions

            use mol_cuda
            implicit none
            integer(4)            :: id, id_int, i, ip
            real(8)               :: prandom
       !     real(8), parameter     :: Half = 0.5d0

            id = (blockidx%x-1)*blockDim%x + threadIDx%x
#if defined (_NORMAL_)
            if (id <= np_d) then
               if (pmove(id) < pspart_d(iptpn_d(id))) then
                  ispartmove_d(id) = 1
                  rotm_d(1,id) = ro_d(1,id) + (ptranx(id)-Half)*dtran_d(iptpn_d(id))
                  rotm_d(2,id) = ro_d(2,id) + (ptrany(id)-Half)*dtran_d(iptpn_d(id))
                  rotm_d(3,id) = ro_d(3,id) + (ptranz(id)-Half)*dtran_d(iptpn_d(id))
                  call PBC_cuda(rotm_d(1,id),rotm_d(2,id),rotm_d(3,id))
               else !if (pmove(id) < pcharge_d(iptpn_d(id))) then
                  ispartmove_d(id) = 2
                  laztm_d(id) = .not. laz_d(id)
                  rotm_d(1,ip) = ro_d(1,ip)
                  rotm_d(2,ip) = ro_d(2,ip)
                  rotm_d(3,ip) = ro_d(3,ip)
                  if (iananweakcharge_d(id) /= 0) laztm_d(iananweakcharge_d(id)) = .not. laz_d(id)
               endif
            end if
#endif
#if defined (_TESTGPU_)
               if (id == 1) then
                  do ip =1, np_d
                     if (pspart_d(iptpn_d(ip)) > One - 1d-14) then
                        ispartmove_d(ip) = 1
                        rotm_d(1,ip) = ro_d(1,ip) + (Random_dev(iseed_trial_d)-Half)*dtran_d(iptpn_d(ip))
                        rotm_d(2,ip) = ro_d(2,ip) + (Random_dev(iseed_trial_d)-Half)*dtran_d(iptpn_d(ip))
                        rotm_d(3,ip) = ro_d(3,ip) + (Random_dev(iseed_trial_d)-Half)*dtran_d(iptpn_d(ip))
                        call PBC_cuda(rotm_d(1,ip),rotm_d(2,ip),rotm_d(3,ip))
                     else
                        prandom = Random_dev(iseed_trial_d)
                        if (prandom < pspart_d(iptpn_d(ip))) then
                           ispartmove_d(ip) = 1
                           rotm_d(1,ip) = ro_d(1,ip) + (Random_dev(iseed_trial_d)-Half)*dtran_d(iptpn_d(ip))
                           rotm_d(2,ip) = ro_d(2,ip) + (Random_dev(iseed_trial_d)-Half)*dtran_d(iptpn_d(ip))
                           rotm_d(3,ip) = ro_d(3,ip) + (Random_dev(iseed_trial_d)-Half)*dtran_d(iptpn_d(ip))
                           call PBC_cuda(rotm_d(1,ip),rotm_d(2,ip),rotm_d(3,ip))
                        else
                           ispartmove_d(ip) = 2
                           laztm_d(ip) = .not. laz_d(ip)
                           rotm_d(1,ip) = ro_d(1,ip)
                           rotm_d(2,ip) = ro_d(2,ip)
                           rotm_d(3,ip) = ro_d(3,ip)
                           if (iananweakcharge_d(ip) /= 0) laztm_d(iananweakcharge_d(ip)) = .not. laz_d(iananweakcharge_d(ip))
                        end if
                     endif
                     !print *, ip, laztm_d(ip), laz_d(ip), ispartmove_d(ip), lcounterion_d(ip)
                  end do
               end if
#endif

      end subroutine CalcNewPositions

      attributes(global) subroutine CalcWeight

         implicit none
         integer(4) :: id

         id = (blockidx%x-1)*blockDim%x + threadIDx%x
         if (ispartmove_d(id) == 1) weight_d(id) = One_d
         if (ispartmove_d(id) == 2) then
            if (laz_d(id)) then
               weight_d(id) = weightd_laz(iptpn_d(id))
            else
               weight_d(id) = weightd_nlaz(iptpn_d(id))
            end if
         endif


      end subroutine CalcWeight


      !! subroutine CalculateUpperPart
      !! running on device, calling from device
      !! Calculates the changes in pair energies which are independent on the acceptance of the trial moves
      !! contains:
      !!  subroutines:
      !!  internal parameters:
      !!                      id: global index of thread and index of particle in global list
      !!  global parameters:
      attributes(global) subroutine CalculateUpperPart(lhsoverlap)

         use precision_m
         implicit none
         integer(4)  :: numblocks
         real(fp_kind),shared :: rox(256)
         real(fp_kind),shared :: roy(256)
         real(fp_kind),shared :: roz(256)
         real(fp_kind),shared :: rojx(256)
         real(fp_kind),shared :: rojy(256)
         real(fp_kind),shared :: rojz(256)
         real(fp_kind),shared :: rotmx(256)
         real(fp_kind),shared :: rotmy(256)
         real(fp_kind),shared :: rotmz(256)
         real(fp_kind) :: rdist
         integer(4) :: iptip_s
         integer(4),shared :: iptjp_s(256)
         real(8)   :: E_s
         real(fp_kind)   :: dx
         real(fp_kind)   :: dy
         real(fp_kind)   :: dz
         real(fp_kind), shared   :: rsumrad_s(16,16)
         integer(4) :: ibuf
         logical, intent(inout)   :: lhsoverlap(np_d)
         integer(4) :: npt_s
         integer(4) :: id, id_int, i, j, jp
         integer(4) :: ictpn_s
         real(fp_kind) :: Etwo_s
         real(fp_kind) :: Etwoold_s
         real(fp_kind) :: Etwonew_s
         real(fp_kind) :: Ebond_s
         real(fp_kind) :: Eclink_s
         real(fp_kind) :: bondk_s
         real(fp_kind) :: bondeq_s
         real(fp_kind) :: bondp_s
         real(fp_kind) :: clinkk_s
         real(fp_kind) :: clinkeq_s
         real(fp_kind) :: clinkp_s
         real(fp_kind) :: usum
         integer(4) :: np_s

               id = ((blockIDx%x-1) * blocksize + threadIDx%x)
               id_int = threadIDx%x
               np_s = np_d
               numblocks = ceiling(real(np_s)/blocksize) -1
               if (id <= np_s) then
                  rotmx(id_int) = rotm_d(1,id)
                  rotmy(id_int) = rotm_d(2,id)
                  rotmz(id_int) = rotm_d(3,id)
                  rox(id_int) = ro_d(1,id)
                  roy(id_int) = ro_d(2,id)
                  roz(id_int) = ro_d(3,id)
                  iptip_s = iptpn_d(id)
                  iptjp_s(id_int) = 0
                  ictpn_s = ictpn_d(id)
                  E_s = 0.0
                  Etwo_s = 0.0
                  Etwoold_s = 0.0
                  Etwonew_s = 0.0
                  Ebond_s = 0.0
                  Eclink_s = 0.0
                  rdist = 0.0
                  lhsoverlap(id) = .false.
                  npt_s = npt_d
                  if ( id_int <= npt_s) then
                     do i=1, npt_s
                        rsumrad_s(id_int,i) = rsumrad(id_int,i)
                     end do
                  end if
               end if


                 call syncthreads

             !! calculate particles that are in other blocks
             do j =blockIDx%x, numblocks
                  call syncthreads
                  if ((id_int + j*blocksize) <= np_s) then
                     rojx(id_int) = ro_d(1,id_int+j*blocksize)
                     rojy(id_int) = ro_d(2,id_int+j*blocksize)
                     rojz(id_int) = ro_d(3,id_int+j*blocksize)
                     iptjp_s(id_int) = iptpn_d(id_int+j*blocksize)
                  end if
                  call syncthreads
                  do i=1, blocksize
                     jp = j*blocksize + i
                     if (jp <= np_s) then
                        !new energy
                        dx = rojx(i) - rotmx(id_int)
                        dy = rojy(i) - rotmy(id_int)
                        dz = rojz(i) - rotmz(id_int)
                        call PBCr2_cuda(dx,dy,dz,rdist)
                        if (rdist < rsumrad_s(iptip_s,iptjp_s(i))) then
                           lhsoverlap(id) = .true.
                        end if
                        !call calcUTabplus(id,jp,rdist,usum)
                        call CallcalcUTabnew(id,jp,rdist,usum)
                        Etwo_d(id) = Etwo_d(id) + usum

                     !old energy
                        dx = rojx(i) - rox(id_int)
                        dy = rojy(i) - roy(id_int)
                        dz = rojz(i) - roz(id_int)
                        call PBCr2_cuda(dx,dy,dz,rdist)
                        call CallcalcUTabold(id,jp,rdist,usum)
                        Etwo_d(id) = Etwo_d(id) - usum
                        Etwoold_d(id) = Etwoold_d(id) + usum
                     end if
                  end do
               end do


              call syncthreads
               !! calculate particles that are in the same block
               do i= 1, blocksize
                     !new energy
                  if (id_int < i) then
                     if ((blockIDx%x-1)*blockDim%x + i <= np_s) then
                        dx = rox(i) - rotmx(id_int)
                        dy = roy(i) - rotmy(id_int)
                        dz = roz(i) - rotmz(id_int)
                        call PBCr2_cuda(dx,dy,dz,rdist)
                        if (rdist < rsumrad_s(iptip_s,iptpn_d(i+(blockIDx%x-1)*blockDim%x))) then
                           lhsoverlap(id) = .true.
                        end if
                        call CallcalcUTabnew(id,i+(blockIDx%x-1)*blockDim%x,rdist,usum)
                        Etwo_d(id) = Etwo_d(id) + usum

                     !old energy
                        dx = rox(i) - rox(id_int)
                        dy = roy(i) - roy(id_int)
                        dz = roz(i) - roz(id_int)
                        call PBCr2_cuda(dx,dy,dz,rdist)
                        call CallcalcUTabold(id,i+(blockIDx%x-1)*blockDim%x,rdist,usum)
                        Etwo_d(id) = Etwo_d(id) - usum
                        Etwoold_d(id) = Etwoold_d(id) + usum
                     end if
                  end if
               end do

               if (id <= np_s) then
                  if (ictpn_s /= 0) then
                        bondk_s = bond_d_k(ictpn_s)
                        bondeq_s = bond_d_eq(ictpn_s)
                        bondp_s = bond_d_p(ictpn_s)
                     do i= 1, 2
                        if (id < bondnn_d(i,id)) then
                           jp = bondnn_d(i,id)
                           dx = ro_d(1, jp) - rotmx(id_int)
                           dy = ro_d(2, jp) - rotmy(id_int)
                           dz = ro_d(3, jp) - rotmz(id_int)
                           call PBCr2_cuda(dx,dy,dz,rdist)
                           Ebond_s = Ebond_s + bondk_s*(sqrt(rdist) - bondeq_s)**bondp_s

                           dx = ro_d(1, jp) - rox(id_int)
                           dy = ro_d(2, jp) - roy(id_int)
                           dz = ro_d(3, jp) - roz(id_int)
                           call PBCr2_cuda(dx,dy,dz,rdist)
                           Ebond_s = Ebond_s - bondk_s*(sqrt(rdist) - bondeq_s)**bondp_s
                        end if
                     end do
                  end if
                  if (lclink_d) then
                     if (nbondcl_d(id) /= 0) then
                           clinkk_s = clink_d_k
                           clinkeq_s = clink_d_eq
                           clinkp_s = clink_d_p
                        do i=1,nbondcl_d(id)
                           if (id < bondcl_d(i,id)) then
                              jp = bondcl_d(i,id)
                              dx = ro_d(1,jp) - rotmx(id_int)
                              dy = ro_d(2,jp) - rotmy(id_int)
                              dz = ro_d(3,jp) - rotmz(id_int)
                              call PBCr2_cuda(dx,dy,dz,rdist)
                              Eclink_s = Eclink_s + clinkk_s*(sqrt(rdist) - clinkeq_s)**clinkp_s

                              dx = ro_d(1,jp) - rox(id_int)
                              dy = ro_d(2,jp) - roy(id_int)
                              dz = ro_d(3,jp) - roz(id_int)
                              call PBCr2_cuda(dx,dy,dz,rdist)
                              Eclink_s = Eclink_s - clinkk_s*(sqrt(rdist) - clinkeq_s)**clinkp_s
                           end if
                        end do
                     end if
                  end if
                  !Etwo_d(id) = Etwo_s
                  !Etwoold_d(id) = Etwoold_s
                  !Etwonew_d(id) = Etwonew_s
                  Ebond_d(id) = Ebond_s
                  Eclink_d(id) = Eclink_s
               end if



      end subroutine CalculateUpperPart





      attributes(grid_global) subroutine MakeDecision_CalcLowerPart(lhsoverlap,ipart,nloop)
         use cooperative_groups
         use precision_m
         use mol_cuda
         implicit none

         real(fp_kind) :: fac_metro
         real(8) :: expmax_d = 87.0d0
         !real(8)  :: numblocks
         real(fp_kind) :: rox
         real(fp_kind) :: roy
         real(fp_kind) :: roz
         real(fp_kind) :: roix
         real(fp_kind) :: roiy
         real(fp_kind) :: roiz
         real(fp_kind) :: rotmx
         real(fp_kind) :: rotmy
         real(fp_kind) :: rotmz
         real(fp_kind) :: rdistold
         real(fp_kind) :: rdistnew
         integer(4) :: iptip_s, iptjp_s
         real(fp_kind)   :: E_s
         real(fp_kind)   :: dx
         real(fp_kind)   :: dy
         real(fp_kind)   :: dz
         real(fp_kind), shared   :: rsumrad_s(16,16)
         logical, intent(inout)   :: lhsoverlap(np_d)
         integer(4), value :: ipart
         integer(4), intent(in) :: nloop
         integer(4) :: npt_s
         integer(4) :: id, id_int, i, j , ibuf, ibuf2, istat
         type(grid_group) :: gg
         integer(4) :: ictpn_s
         real(fp_kind) :: Etwo_s
         real(fp_kind) :: Etwoold_s
         real(fp_kind) :: Etwonew_s
         real(fp_kind) :: Ebond_s
         real(fp_kind) :: Eclink_s
         real(fp_kind) :: bondk_s
         real(fp_kind) :: bondeq_s
         real(fp_kind) :: bondp_s
         real(fp_kind) :: clinkk_s
         real(fp_kind) :: clinkeq_s
         real(fp_kind) :: clinkp_s
         real(fp_kind) :: dured
         real(fp_kind) :: usum, d
         integer(4), ipartmin, ipartmax, idecision
         integer(4), np_s
         real(fp_kind) :: beta_s
         real(fp_kind), shared :: du_s
         real(fp_kind), shared :: dutwo_s
         real(fp_kind), shared :: dubond_s
         real(fp_kind), shared :: duclink_s


               gg = this_grid()
               np_s = np_d
               ipartmin = blocksize*iblock2_d*(ipart-1)!np_s/nloop*(ipart-1)
               if (ipart == nloop) then
                  ipartmax = np_s
               else
                  ipartmax = iblock2_d*blocksize*ipart!np_s/nloop*ipart
               end if
               id = ((blockIDx%x-1) * blocksize + threadIDx%x)+ipartmin
               id_int = threadIDx%x

               if (id <= np_s) then
                  rotmx = rotm_d(1,id)
                  rotmy = rotm_d(2,id)
                  rotmz = rotm_d(3,id)
                  rox = ro_d(1,id)
                  roy = ro_d(2,id)
                  roz = ro_d(3,id)
                  iptip_s = iptpn_d(id)
                  E_s = 0.0
                  !Etwo_s = Etwo_d(id)
                  !Etwoold_s = Etwoold_d(id)
                  !Etwonew_s = Etwonew_d(id)
                  Ebond_s = Ebond_d(id)
                  Eclink_s = Eclink_d(id)
                    !numblocks = np_s / blocksize
                  npt_s = npt_d
                  du_s = 0.0
                  dutwo_s = 0.0
                  dubond_s = 0.0
                  duclink_s = 0.0
                  beta_s = beta_d
                  ictpn_s = ictpn_d(id)
                  if ( id_int <= npt_s) then
                     do i=1, npt_s
                        rsumrad_s(id_int,i) = rsumrad(id_int,i)
                     end do
                  end if
                  if (ictpn_s /= 0) then
                           bondk_s = bond_d_k(ictpn_s)
                           bondeq_s = bond_d_eq(ictpn_s)
                           bondp_s = bond_d_p(ictpn_s)
                  end if
                  if (nbondcl_d(id) /= 0) then
                           clinkk_s = clink_d_k
                           clinkeq_s = clink_d_eq
                           clinkp_s = clink_d_p
                  end if
               end if
               call syncthreads

         do i = 1+ipartmin, ipartmax
            !if (id == iananweakcharge_d(i)) Etwoold_d(id) = Etwoold_s
            !call syncthreads(gg)
            if (id == i) then
               if (ispartmove_d(id) == 2) then
                  if ( laztm_d(id)) then
                     !Etwo_s = Etwo_s + Etwoold_d(iananweakcharge_d(id))
                     Etwo_d(id) = Etwo_d(id) + Etwoold_d(iananweakcharge_d(id))
                  else
                     !Etwo_s = Etwo_s - Etwoold_d(iananweakcharge_d(id))
                     Etwo_d(id) = Etwo_d(id) - Etwoold_d(iananweakcharge_d(id))
                  end if
               end if
               if (lcounterion_d(id)) then
                  if (laz_d(id) == .false.) Etwo_d(id) = 0.0
               end if

                  E_s = Etwo_d(id) + Ebond_s + Eclink_s
                  dured = beta_s*E_s
                  if (lhsoverlap(id) == .true.) then
                      idecision = 4   !imchsreject
                  else
                     if (dured > expmax_d) then
                        idecision = 2 ! energy rejected
                     else if (dured < -expmax_d) then
                        idecision = 1   !accepted
                     else
                        fac_metro = weight_d(id)*exp(-dured)
                        if (fac_metro > One_d) then
                           idecision = 1 !accepted
#if defined (_TESTGPU_)
                        else if (fac_metro > Random_dev2(iseed2_d)) then
#else
                        else if (fac_metro > pmetro(id)) then
#endif
                           idecision = 1 ! accepted
                        else
                           idecision = 2 ! energy rejected
                        end if
                     end if
                     if (idecision == 1) then
                        if (ispartmove_d(id) == 1) then
                           ro_d(1,id) = rotmx
                           ro_d(2,id) = rotmy
                           ro_d(3,id) = rotmz
                        else
                           if (ispartmove_d(id) == 2) then
                              if (lweakcharge_d) then
                                 laz_d(id) = laztm_d(id)
                                 if (iananweakcharge_d(id) /= 0) then
                                    laz_d(iananweakcharge_d(id)) = laztm_d(iananweakcharge_d(id))
                                 end if
                              end if
                           end if
                        end if
                           dutot_d = dutot_d + E_s
                           dutwo_d = dutwo_d + Etwo_d(id)
                           dubond_d = dubond_d + Ebond_s
                           duclink_d = duclink_d + Eclink_s
                     end if
                  end if
                  mcstat_d(iptip_s,idecision) = mcstat_d(iptip_s,idecision) + 1
            end if
            call syncthreads(gg)
            !if (lweakcharge_d) then
            !   if (id == iananweakcharge_d(i)) then
            !      if (laz_d(id) == .false.) Etwo_s = 0.0
            !   end if
            !end if
            roix = ro_d(1,i)
            roiy = ro_d(2,i)
            roiz = ro_d(3,i)
            iptjp_s = iptpn_d(i)
         if ( id <= np_s) then
            if ( id > i) then
                     dx = roix - rotmx
                     dy = roiy - rotmy
                     dz = roiz - rotmz
                     call PBCr2_cuda(dx,dy,dz,rdistnew)
                     if (rdistnew < rsumrad_s(iptip_s,iptjp_s)) then
                        lhsoverlap(id) = .true.
                     end if
                     call CallcalcUTabnew(id,i,rdistnew,usum)
                     Etwo_d(id) = Etwo_d(id) + usum
                  !old energy
                     dx = roix - rox
                     dy = roiy - roy
                     dz = roiz - roz
                     call PBCr2_cuda(dx,dy,dz,rdistold)
                     call CallcalcUTabold(id,i,rdistold,usum)
                     Etwo_d(id) = Etwo_d(id) - usum
                     Etwoold_d(id) = Etwoold_d(id) + usum

               !calculate bonds
               if (ictpn_d(i) /= 0) then
                  do j= 1, 2
                     if (i == bondnn_d(j,id)) then
                        Ebond_s = Ebond_s + bondk_s*((sqrt(rdistnew)-bondeq_s)**bondp_s - &
                           (sqrt(rdistold) - bondeq_s)**bondp_s)
                     end if

                  end do
               end if
               ! calculate crosslinks
               if (lclink_d) then
                  do j=1,nbondcl_d(i)
                     if (id == bondcl_d(j,i)) then
                        Ebond_s = Ebond_s + clinkk_s*((sqrt(rdistnew)-clinkeq_s)**clinkp_s - &
                           (sqrt(rdistold) - clinkeq_s)**clinkp_s)
                     end if
                  end do
               end if
            end if
         end if
            call syncthreads(gg)
      end do


      end subroutine MakeDecision_CalcLowerPart


      attributes(global) subroutine CalcLowerPart2(lhsoverlap, ipart, nloop)

         !use cooperative_groups
         use precision_m
         implicit none
         integer(4)  ::numblocks
         integer(4)  ::numblocks_old
         real(fp_kind),shared :: rox(256)
         real(fp_kind),shared :: roy(256)
         real(fp_kind),shared :: roz(256)
         real(fp_kind),shared :: rojx(256)
         real(fp_kind),shared :: rojy(256)
         real(fp_kind),shared :: rojz(256)
         real(fp_kind),shared :: rotmx(256)
         real(fp_kind),shared :: rotmy(256)
         real(fp_kind),shared :: rotmz(256)
         real(fp_kind) :: rdist
         integer(4) :: iptip_s
         integer(4) :: iptjp_s(256)
         real(fp_kind)   :: E_s
         real(fp_kind) :: Etwo_s
         real(fp_kind) :: Etwoold_s
         real(fp_kind) :: Etwonew_s
         real(fp_kind) :: Ebond_s
         real(fp_kind) :: Eclink_s
         real(fp_kind)   :: dx
         real(fp_kind)   :: dy
         real(fp_kind)   :: dz
         real(fp_kind), shared   :: rsumrad_s(16,16)
         logical, intent(inout)   :: lhsoverlap(np_d)
         integer(4) :: npt_s
         integer(4) :: id, id_int!,
         integer(4) :: i, j, q, jp!, ibuf
         integer(4), value :: ipart
         integer(4), intent(in) :: nloop
         integer(4) :: ipart_2, np_s
         integer(4) :: ictpn_s
         real(fp_kind) :: bondk_s
         real(fp_kind) :: bondeq_s
         real(fp_kind) :: bondp_s
         real(fp_kind) :: clinkk_s
         real(fp_kind) :: clinkeq_s
         real(fp_kind) :: clinkp_s
         real(fp_kind) :: usum

         np_s = np_d
         numblocks = iblock2_d * ipart
         numblocks_old = iblock2_d * (ipart - 1) + 1
         npt_s = npt_d

         do ipart_2 = ipart, nloop - 1
            id = ((blockIDx%x-1) * blockDim%x + threadIDx%x)+(iblock2_d*blockDim%x*ipart_2)
               id_int = threadIDx%x
               if (id <= np_s) then
                  rotmx(id_int) = rotm_d(1,id)
                  rotmy(id_int) = rotm_d(2,id)
                  rotmz(id_int) = rotm_d(3,id)
                  rox(id_int) = ro_d(1,id)
                  roy(id_int) = ro_d(2,id)
                  roz(id_int) = ro_d(3,id)
                  iptip_s = iptpn_d(id)
                  iptjp_s(id_int) = 0
                  !E_s = E_g(id)
                  !Etwo_s = Etwo_d(id)
                  !Etwoold_s = Etwoold_d(id)
                  !Etwonew_s = Etwonew_d(id)
                  Ebond_s = Ebond_d(id)
                  Eclink_s = Eclink_d(id)
                  iptip_s = iptpn_d(id)
                  ictpn_s = ictpn_d(id)
                  do j = 1, npt_s
                     do i=1, npt_s
                        rsumrad_s(j,i) = rsumrad(j,i)
                     end do
                  end do
                  if (ictpn_s /= 0) then
                        bondk_s = bond_d_k(ictpn_s)
                        bondeq_s = bond_d_eq(ictpn_s)
                        bondp_s = bond_d_p(ictpn_s)
                  end if
               end if
               call syncthreads




            do j = numblocks_old, numblocks
               call syncthreads
                  rojx(id_int) = ro_d(1,id_int+(j-1)*blockDim%x)
                  rojy(id_int) = ro_d(2,id_int+(j-1)*blockDim%x)
                  rojz(id_int) = ro_d(3,id_int+(j-1)*blockDim%x)
                  iptjp_s(id_int) = iptpn_d(id_int+(j-1)*blockDim%x)
               call syncthreads
               if (id <= np_s) then
                  do i=1, blocksize
                     !new energy
                        dx = rojx(i) - rotmx(id_int)
                        dy = rojy(i) - rotmy(id_int)
                        dz = rojz(i) - rotmz(id_int)
                        !dx = ro_d(i+(j-1)*blocksize) - ro_d(id)
                        !dy = ro_d(i+(j-1)*blocksize) - ro_d(id)
                        !dz = ro_d(i+(j-1)*blocksize) - ro_d(id)
                        call PBCr2_cuda(dx,dy,dz,rdist)
                        !if (rdist < rsumrad_s(iptip_s,iptjp_s(i))) then
                        if (rdist < rsumrad_s(iptip_s,iptpn_d(i+(j-1)*blockDim%x))) then
                           lhsoverlap(id) = .true.
                        end if
                        !call calcUTabplus(id,(j-1)*blockDim%x + i,rdist,usum)
                        call CallcalcUTabnew(id,(j-1)*blockDim%x + i,rdist,usum)
                        !E_s = E_s + usum
                        !Etwo_s = Etwo_s + usum
                        Etwo_d(id) = Etwo_d(id) + usum
                        !Etwonew_s = Etwonew_s + usum


                     !old energy
                        dx = rojx(i) - rox(id_int)
                        dy = rojy(i) - roy(id_int)
                        dz = rojz(i) - roz(id_int)
                        call PBCr2_cuda(dx,dy,dz,rdist)
                        !call calcUTabminus(id,(j-1)*blockDim%x + i,rdist,usum)
                        call CallcalcUTabold(id,(j-1)*blockDim%x + i,rdist,usum)
                        !E_s = E_s - usum
                        !Etwo_s = Etwo_s - usum
                        Etwo_d(id) = Etwo_d(id) - usum
                        !Etwoold_s = Etwoold_s + usum
                        Etwoold_d(id) = Etwoold_d(id) + usum
                  end do
               end if
            end do
            if (id <= np_s) then
               if (ictpn_s /= 0) then
                  do q= 1, 2
                     if ((numblocks_old-1)*blockDim%x < bondnn_d(q,id) .and. numblocks*blockDim%x >= bondnn_d(q,id) ) then
                        jp = bondnn_d(q,id)
                        dx = ro_d(1,jp) - rotmx(id_int)
                        dy = ro_d(2,jp) - rotmy(id_int)
                        dz = ro_d(3,jp) - rotmz(id_int)
                        call PBCr2_cuda(dx,dy,dz,rdist)
                        rdist = sqrt(rdist)
                        !E_s = E_s + bondk_s*(rdist - bondeq_s)**bondp_s
                        Ebond_s = Ebond_s + bondk_s*(rdist - bondeq_s)**bondp_s

                        dx = ro_d(1,jp) - rox(id_int)
                        dy = ro_d(2,jp) - roy(id_int)
                        dz = ro_d(3,jp) - roz(id_int)
                        call PBCr2_cuda(dx,dy,dz,rdist)
                        rdist = sqrt(rdist)
                        !E_s = E_s - bondk_s*(rdist - bondeq_s)**bondp_s
                        Ebond_s = Ebond_s - bondk_s*(rdist - bondeq_s)**bondp_s
                     end if
                  end do
               end if
               if (lclink_d) then
                  if (nbondcl_d(id) /= 0) then
                        clinkk_s = clink_d_k
                        clinkeq_s = clink_d_eq
                        clinkp_s = clink_d_p
                     do i=1,nbondcl_d(id)
                        if ((numblocks_old-1)*blockDim%x < bondcl_d(i,id).and. (numblocks*blockDim%x >= bondcl_d(i,id))) then
                           jp = bondcl_d(i,id)
                           dx = ro_d(1,jp) - rotmx(id_int)
                           dy = ro_d(2,jp) - rotmy(id_int)
                           dz = ro_d(3,jp) - rotmz(id_int)
                           call PBCr2_cuda(dx,dy,dz,rdist)
                           !E_s = E_s + clinkk_s*(sqrt(rdist) - clinkeq_s)**clinkp_s
                           Eclink_s = Eclink_s + clinkk_s*(sqrt(rdist) - clinkeq_s)**clinkp_s

                           dx = ro_d(1,jp) - rox(id_int)
                           dy = ro_d(2,jp) - roy(id_int)
                           dz = ro_d(3,jp) - roz(id_int)
                           call PBCr2_cuda(dx,dy,dz,rdist)
                           !E_s = E_s - clinkk_s*(sqrt(rdist) - clinkeq_s)**clinkp_s
                           Eclink_s = Eclink_s - clinkk_s*(sqrt(rdist) - clinkeq_s)**clinkp_s
                        end if
                     end do
                  end if
               end if
               !E_g(id) = E_s
               !Etwo_d(id) = Etwo_s
               !Etwoold_d(id) = Etwoold_s
               !Etwonew_d(id) = Etwonew_s
               Ebond_d(id) = Ebond_s
               Eclink_d(id) = Eclink_s
            end if
               call syncthreads
         end do
      end subroutine CalcLowerPart2



      attributes(device) subroutine CallcalcUTabnew(id,i,rdist,E)

      implicit none
      integer(4), intent(in) :: id
      integer(4), intent(in) :: i
      real(fp_kind), intent(inout) :: E
      real(fp_kind), intent(in) :: rdist

      if(lcharge_d) call calcUTab(id,i,rdist,E)
      if(lweakcharge_d) call calcUTabnew_weak(id,i,rdist,E)


      end subroutine CallcalcUTabnew

      attributes(device) subroutine CallcalcUTabold(id,i,rdist,E)

      implicit none
      integer(4), intent(in) :: id
      integer(4), intent(in) :: i
      real(fp_kind), intent(inout) :: E
      real(fp_kind), intent(in) :: rdist

      if(lcharge_d) call calcUTab(id,i,rdist,E)
      if(lweakcharge_d) call calcUTabold_weak(id,i,rdist,E)


      end subroutine CallcalcUTabold

      attributes(device) subroutine calcUTab(id,i,rdist,E)

      implicit none
      integer(4), intent(in) :: id
      integer(4), intent(in) :: i
      real(fp_kind), intent(inout) :: E
      real(fp_kind), intent(in) :: rdist
      integer(4) :: ibuf
      real(fp_kind) :: d

        ibuf = iubuflow_d(iptpt_d(iptpn_d(id),iptpn_d(i)))
        !ibuf = 1
         do
            if (rdist >= ubuf_d(ibuf)) exit
            ibuf = ibuf+12
         end do
            d = rdist - ubuf_d(ibuf)
         E = ubuf_d(ibuf+1)+d*(ubuf_d(ibuf+2)+d*(ubuf_d(ibuf+3)+ &
               d*(ubuf_d(ibuf+4)+d*(ubuf_d(ibuf+5)+d*ubuf_d(ibuf+6)))))




      end subroutine calcUTab

      attributes(device) subroutine calcUTabnew_weak(id,i,rdist,E)

      implicit none
      integer(4), intent(in) :: id
      integer(4), intent(in) :: i
      real(fp_kind), intent(inout) :: E
      real(fp_kind), intent(in) :: rdist
      logical :: lazj
      integer(4) :: ibuf
      real(fp_kind) :: d

      !print *, id, laztm_d(id), i, laz_d(i)
      if (laztm_d(id) .OR. lcounterion_d(id)) then
         if (i == iananweakcharge_d(id)) then
            lazj = laztm_d(i)
         else
            lazj = laz_d(i)
         end if
      !print *, id, laztm_d(id), lazj, i
         if (lazj == .true.) then
           ibuf = iubuflow_d(iptpt_d(iptpn_d(id),iptpn_d(i)))
            do
               if (rdist >= ubuf_d(ibuf)) exit
               ibuf = ibuf+12
            end do
               d = rdist - ubuf_d(ibuf)
            E = ubuf_d(ibuf+1)+d*(ubuf_d(ibuf+2)+d*(ubuf_d(ibuf+3)+ &
                  d*(ubuf_d(ibuf+4)+d*(ubuf_d(ibuf+5)+d*ubuf_d(ibuf+6)))))
         else
               E = Zero_d
         end if
      else
         E = Zero_d
      end if
      !print *, id, i, E

      end subroutine calcUTabnew_weak

      attributes(device) subroutine calcUTabold_weak(id,i,rdist,E)

      implicit none
      integer(4), intent(in) :: id
      integer(4), intent(in) :: i
      real(fp_kind), intent(inout) :: E
      real(fp_kind), intent(in) :: rdist
      logical  :: lcharge
      integer(4) :: ibuf
      real(fp_kind) :: d

      if (laz_d(id) .OR. lcounterion_d(id)) then
         if (laz_d(i)) then
           ibuf = iubuflow_d(iptpt_d(iptpn_d(id),iptpn_d(i)))
            do
               if (rdist >= ubuf_d(ibuf)) exit
               ibuf = ibuf+12
            end do
               d = rdist - ubuf_d(ibuf)
            E = ubuf_d(ibuf+1)+d*(ubuf_d(ibuf+2)+d*(ubuf_d(ibuf+3)+ &
                  d*(ubuf_d(ibuf+4)+d*(ubuf_d(ibuf+5)+d*ubuf_d(ibuf+6)))))
         else
            E = Zero_d
         end if
      else
         E = Zero_d
      end if

      end subroutine calcUTabold_weak


subroutine startUTwoBodyAAll(lhsoverlap)

   use mol_cuda
   use cudafor
   use precision_m

   implicit none
   logical, intent(inout) :: lhsoverlap
   integer(4) :: numblocks
   integer(4) :: sizeofblocks =512
   integer(4) :: isharedmem


           call TransferDUTotalVarToDevice
           numblocks = floor(Real((nptm*np)/sizeofblocks)) + 1
           lhsoverlap = .false.
           isharedmem = 2*sizeofblocks*fp_kind + sizeofblocks*4 + threadssum*(nptpt+1)*fp_kind
           lhsoverlap_d = .false.

           call UTwoBodyAAll<<<numblocks,sizeofblocks,isharedmem>>>(lhsoverlap_d)                ! calculate new two-body potential energy
           lhsoverlap = lhsoverlap_d
end subroutine startUTwoBodyAAll

attributes(global) subroutine UTwoBodyAAll(lhsoverlap)


   !use EnergyModule
   use mol_cuda
   use cudafor
   use precision_m
   implicit none

   logical,    intent(out) :: lhsoverlap
   !character(40), parameter :: txroutine ='UTwoBodyANew'

   integer(4) :: ip, iploc, ipt, jploc, jpt, iptjpt, ibuf,jp, i, j
   real(fp_kind)    :: dx, dy, dz, r2, d
   integer(4) :: tidx, t, tidx_int, istat
   integer(4),shared :: iptjpt_arr(blockDim%x)
   real(fp_kind), shared :: usum1(blockDim%x), usum2(blockDim%x)
   real(fp_kind), shared ::  usum_aux1(threadssum_d,0:nptpt_d)
!   logical    :: EllipsoidOverlap, SuperballOverlap
   tidx = blockDim%x * (blockIdx%x - 1) + threadIdx%x  !global thread index 1 ...
   tidx_int = threadIDx%x
   iploc = ceiling(real((tidx-1)/np_d))+1
   if (iploc <= nptm_d) then
      ip = ipnptm_d(iploc)
   end if
   jp = mod(tidx-1,np_d)+1
     usum_aux1 = 0.0
    iptjpt = 0
    iptjpt_arr(tidx_int) = iptjpt
    usum1(tidx_int) = 0.0
   call syncthreads


   if (tidx <= nptm_d*np_d) then
      ipt = iptpn_d(ip)
        if ( ip /= jp ) then
             jpt = iptpn_d(jp)
             iptjpt = iptpt_d(ipt,jpt)
             iptjpt_arr(tidx_int) = iptjpt
             if (.not. lptm_d(jp)) then
                dx = rotm_d(1,iploc)-ro_d(1,jp)
                dy = rotm_d(2,iploc)-ro_d(2,jp)
                dz = rotm_d(3,iploc)-ro_d(3,jp)
             else
                if (ip < jp) then
                  do jploc = 1, nptm_d
                     dx = rotm_d(1,iploc)-rotm_d(1,jploc)
                     dy = rotm_d(2,iploc)-rotm_d(2,jploc)
                     dz = rotm_d(3,iploc)-rotm_d(3,jploc)
                   end do
                 else
                      goto 400
                 end if
              end if
              call PBCr2_cuda(dx,dy,dz,r2)
              if (lellipsoid_d) Then
              ! if (EllipsoidOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
              end if
              if (lsuperball_d) Then
              ! if (SuperballOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp))) goto 400
              end if
              if (r2 > rcut2_d) then
                !do not anything
              else if (r2 < r2atat_d(iptjpt))then
                  lhsoverlap = .true.
              else if (r2 < r2umin_d(iptjpt))then       ! outside lower end
                  lhsoverlap = .true.
              else
                 ibuf = iubuflow_d(iptjpt)
                 do
                    if (r2 >= ubuf_d(ibuf)) exit
                       ibuf = ibuf+12
                    end do
                 d = r2-ubuf_d(ibuf)
                 usum1(tidx_int) = ubuf_d(ibuf+1)+d*(ubuf_d(ibuf+2)+d*(ubuf_d(ibuf+3)+ &
                              d*(ubuf_d(ibuf+4)+d*(ubuf_d(ibuf+5)+d*ubuf_d(ibuf+6)))))
              end if
              if (.not. lptm_d(jp)) then
               dx = ro_d(1,ip)-ro_d(1,jp)
               dy = ro_d(2,ip)-ro_d(2,jp)
               dz = ro_d(3,ip)-ro_d(3,jp)
              else if (ip < jp) then
                do jploc = 1, nptm_d
                  dx = ro_d(1,ip)-ro_d(1,jp)
                  dy = ro_d(2,ip)-ro_d(2,jp)
                  dz = ro_d(3,ip)-ro_d(3,jp)
                end do
              else
                  goto 400
              end if
               call PBCr2_cuda(dx,dy,dz,r2)
               if (lellipsoid_d) Then
                 ! if (EllipsoidOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
               end if
               if (lsuperball_d) Then
                 ! if (SuperballOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp))) goto 400
               end if
               if (r2 > rcut2_d) goto 400
                   !do not anything
               if (r2 < r2atat_d(iptjpt)) goto 400
                    ! lhsoverlap = .true.
               if (r2 < r2umin_d(iptjpt)) goto 400      ! outside lower end
                    ! lhsoverlap = .true.
              ibuf = iubuflow_d(iptjpt)
              do
                 if (r2 >= ubuf_d(ibuf)) exit
                 ibuf = ibuf+12
              end do
              d = r2-ubuf_d(ibuf)
              usum2(tidx_int) = ubuf_d(ibuf+1)+d*(ubuf_d(ibuf+2)+d*(ubuf_d(ibuf+3)+ &
                           d*(ubuf_d(ibuf+4)+d*(ubuf_d(ibuf+5)+d*ubuf_d(ibuf+6)))))
                        usum1(tidx_int) = usum1(tidx_int) - usum2(tidx_int)
        end if
     end if

  400 continue
       call syncthreads
       if (tidx_int <= threadssum_d) then
          do i = 1, blockDim%x/threadssum_d
             usum_aux1(tidx_int,iptjpt_arr(threadssum_d*(i-1) + tidx_int)) = &
                usum_aux1(tidx_int,iptjpt_arr(threadssum_d*(i-1) + tidx_int)) + usum1(threadssum_d*(i-1)+tidx_int)
          end do
       end if
          call syncthreads
       if (tidx_int == 1) then
          do i = 2, threadssum_d
             do j = 1, nptpt_d
            usum_aux1(1,j) = usum_aux1(1,j) + usum_aux1(i,j)
            end do
          end do

          do i = 1, nptpt_d
             istat = atomicAdd(utwobnew_d(i),usum_aux1(1,i))
          end do
       end if

end subroutine UTwoBodyAAll

attributes(grid_global) subroutine MCPass_cuda

   use mol_cuda
   use cooperative_groups
   !use EnergyModule
   use cudafor
   implicit none
   integer(4) :: np_s
   type(grid_group) :: gg
   logical :: lhsoverlap
   integer(4) :: ip, iploc, ipt, jploc, jpt, iptjpt, ibuf,jp, i, j, n,m, idecision,iloops_s,ipartloop
   real(fp_kind) :: beta_s, fac_metro, dured
   real(fp_kind)    :: dx, dy, dz, r2new, r2old, d
   integer(4) :: tidx, t, tidx_int, istat
   integer(4),shared :: iptjpt_arr(blockDim%x)
   real(8) :: expmax_d = 87.0d0
   real(fp_kind), shared :: usum1(blockDim%x), usum2(blockDim%x)
   real(fp_kind), shared ::  usum_aux1(threadssum_d,0:nptpt_d)

   gg = this_grid()
   np_s = np_d
   beta_s = beta_d
   iloops_s = iloops_d
   ipartloop = blockDim%x*iblock2_d

   do n = 1, np_s
       call syncthreads(gg)
       lhsoverlap_d = .false.
       dubond_d = 0.0
       duclink_d = 0.0
       utwobnew_d = 0.0
       usum_aux1 = 0.0
       call syncthreads(gg)
      do m = 1, iloops_s
         tidx = (blockDim%x * (blockIdx%x - 1)) + ipartloop*(m-1) + threadIdx%x  !global thread index 1 ...
         tidx_int = threadIDx%x
         !iploc = ceiling(real((tidx-1)/np_s))+1
         !if (iploc <= nptm_d) then
         !   ip = ipnptm_d(iploc)
         !end if
         !jp = mod(tidx-1,np_s)+1
         jp = tidx
         iptjpt = 0
         iptjpt_arr(tidx_int) = iptjpt
         usum1(tidx_int) = 0.0
         call syncthreads


         if (tidx <= np_s) then
            ipt = iptpn_d(n)
              if ( jp /= n ) then
                   jpt = iptpn_d(jp)
                   iptjpt = iptpt_d(ipt,jpt)
                   iptjpt_arr(tidx_int) = iptjpt
                      dx = rotm_d(1,n)-ro_d(1,jp)
                      dy = rotm_d(2,n)-ro_d(2,jp)
                      dz = rotm_d(3,n)-ro_d(3,jp)
                    call PBCr2_cuda(dx,dy,dz,r2new)
                    !if (lellipsoid_d) Then
                    ! if (EllipsoidOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
                    !end if
                    !if (lsuperball_d) Then
                    ! if (SuperballOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp))) goto 400
                    !end if
                    if (r2new > rcut2_d) then
                      !do not anything
                    else if (r2new < r2atat_d(iptjpt))then
                        lhsoverlap_d = .true.
                    else if (r2new < r2umin_d(iptjpt))then       ! outside lower end
                        lhsoverlap_d = .true.
                    else
                       ibuf = iubuflow_d(iptjpt)
                       do
                          if (r2new >= ubuf_d(ibuf)) exit
                             ibuf = ibuf+12
                       end do
                       d = r2new-ubuf_d(ibuf)
                       usum1(tidx_int) = ubuf_d(ibuf+1)+d*(ubuf_d(ibuf+2)+d*(ubuf_d(ibuf+3)+ &
                                    d*(ubuf_d(ibuf+4)+d*(ubuf_d(ibuf+5)+d*ubuf_d(ibuf+6)))))
                    end if

                    !old energy
                     dx = ro_d(1,n)-ro_d(1,jp)
                     dy = ro_d(2,n)-ro_d(2,jp)
                     dz = ro_d(3,n)-ro_d(3,jp)
                     call PBCr2_cuda(dx,dy,dz,r2old)
                     !if (lellipsoid_d) Then
                       ! if (EllipsoidOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
                     !end if
                     !if (lsuperball_d) Then
                       ! if (SuperballOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp))) goto 400
                     !end if
                     if (r2old > rcut2_d) goto 400
                         !do not anything
                     if (r2old < r2atat_d(iptjpt)) goto 400
                          ! lhsoverlap = .true.
                     if (r2old < r2umin_d(iptjpt)) goto 400      ! outside lower end
                          ! lhsoverlap = .true.
                    ibuf = iubuflow_d(iptjpt)
                    do
                       if (r2old >= ubuf_d(ibuf)) exit
                       ibuf = ibuf+12
                    end do
                    d = r2old-ubuf_d(ibuf)
                    usum2(tidx_int) = ubuf_d(ibuf+1)+d*(ubuf_d(ibuf+2)+d*(ubuf_d(ibuf+3)+ &
                                 d*(ubuf_d(ibuf+4)+d*(ubuf_d(ibuf+5)+d*ubuf_d(ibuf+6)))))
                    usum1(tidx_int) = usum1(tidx_int) - usum2(tidx_int)
              end if
           end if

        400 continue
             call syncthreads
             if (tidx_int <= threadssum_d) then
                do i = 1, blockDim%x/threadssum_d
                   usum_aux1(tidx_int,iptjpt_arr(threadssum_d*(i-1) + tidx_int)) = &
                      usum_aux1(tidx_int,iptjpt_arr(threadssum_d*(i-1) + tidx_int)) + usum1(threadssum_d*(i-1)+tidx_int)
                end do
             end if
             call syncthreads
                  !calculate bonds
              if (tidx <= np_s) then
                  if (ictpn_d(n) /= 0) then
                     do j= 1, 2
                        if (n == bondnn_d(j,tidx)) then

                           usum1(tidx_int) = bond_d_k(ictpn_d(n))*((sqrt(r2new)-bond_d_eq(ictpn_d(n)))**bond_d_p(ictpn_d(n)) - &
                                       (sqrt(r2old) - bond_d_eq(ictpn_d(n)))**bond_d_p(ictpn_d(n)))
                           istat = atomicAdd(dubond_d,usum1(tidx_int))
                        end if
                     end do
                  end if
                  ! calculate crosslinks
                  if (lclink_d) then
                     do j=1,nbondcl_d(n)
                        if (tidx == bondcl_d(j,n)) then

                           usum1(tidx_int) = clink_d_k*((sqrt(r2new)-clink_d_eq)**clink_d_p - &
                              (sqrt(r2old) - clink_d_eq)**clink_d_p)
                           istat = atomicAdd(duclink_d,usum1(tidx_int))
                        end if
                     end do
                  end if
               end if
         end do
         !reduction of dutwobody
          call syncthreads
          if (tidx_int == 1) then
             do i = 2, threadssum_d
                do j = 1, nptpt_d
               usum_aux1(1,j) = usum_aux1(1,j) + usum_aux1(i,j)
               end do
             end do

             do i = 1, nptpt_d
                istat = atomicAdd(utwobnew_d(i),usum_aux1(1,i))
             end do
          end if
         ! call syncthreads(gg)
          !if (tidx == np_s) then
          !   do i = 1, nptpt_d
          !      istat = atomicAdd(utwobnew_d(0), utwobnew_d(i))
          !   end do
          !end if

              call syncthreads(gg)
               if (tidx == np_s) then
                     do i = 1, nptpt_d
                        istat = atomicAdd(utwobnew_d(0), utwobnew_d(i))
                     end do
                     dutot_d = utwobnew_d(0) + dubond_d + duclink_d
                     dured = beta_s*dutot_d
                     if (lhsoverlap_d == .true.) then
                         idecision = 4   !imchsreject
                     else
                        if (dured > expmax_d) then
                           idecision = 2 ! energy rejected
                        else if (dured < -expmax_d) then
                           idecision = 1   !accepted
                        else
                           fac_metro = exp(-dured)
                           if (fac_metro > One_d) then
                              idecision = 1 !accepted
                           else if (fac_metro > Random_dev2(iseed2_d)) then
                           !else if (fac_metro > pmetro(id)) then
                              idecision = 1 ! accepted
                           else
                              idecision = 2 ! energy rejected
                           end if
                        end if
                        if (idecision == 1) then
                              ro_d(1,n) = rotm_d(1,n)
                              ro_d(2,n) = rotm_d(2,n)
                              ro_d(3,n) = rotm_d(3,n)
                              !utot_d = utot_d + E_s
                              utot_d = utot_d + dutot_d
                        end if
                     end if
                     mcstat_d(iptpn_d(n),idecision) = mcstat_d(iptpn_d(n),idecision) + 1
               end if
              call syncthreads(gg)
   end do

end subroutine MCPass_cuda


!************************************************************************
!> \page PBCr2_cuda
!! **PBCr2_cuda**
!! *apply periodic boundary conditions and calculate r**2 on GPU*
!************************************************************************


attributes(device) subroutine PBCr2_cuda(dx,dy,dz,r2)!,boxlen,boxlen2,lPBC,lbcbox,&
                              !lbcrd,lbcto)

   implicit none

   real(fp_kind), intent(inout) :: dx, dy, dz
   real(fp_kind), intent(inout) :: r2
!   real(8), intent(inout)  :: boxlen, boxlen2
!   logical, intent(inout)  :: lPBC, lbcbox, lbcrd, lbcto
   !real(8)              :: Threehalf = 1.5d0
   !real(8)              :: SqTwo = sqrt(2.0d0)

   if (lPBC_d) then                                                              ! periodic boundary condition
      if (lbcbox_d) then                                                         ! box-like cell
         if (abs(dx) > boxlen2_d(1)) dx = dx - sign(dpbc_d(1),dx)
         if (abs(dy) > boxlen2_d(2)) dy = dy - sign(dpbc_d(2),dy)
         if (abs(dz) > boxlen2_d(3)) dz = dz - sign(dpbc_d(3),dz)
      else if (lbcrd_d) then                                                     ! rhombic dodecahedral cell
         if (abs(dx) > boxlen2_d(1)) dx = dx - sign(boxlen_d(1),dx)
         if (abs(dy) > boxlen2_d(2)) dy = dy - sign(boxlen_d(2),dy)
         if (abs(dz) > boxlen2_d(3)) dz = dz - sign(boxlen_d(3),dz)
         if (abs(dx) + abs(dy) + SqTwo_d*abs(dz) > boxlen_d(1)) then
            dx = dx - sign(boxlen2_d(1),dx)
            dy = dy - sign(boxlen2_d(2),dy)
            dz = dz - sign(boxlen2_d(3),dz)
         end if
      else if (lbcto_d) then                                                     ! truncated octahedral cell
         if (abs(dx) > boxlen2_d(1)) dx = dx - sign(boxlen_d(1),dx)
         if (abs(dy) > boxlen2_d(2)) dy = dy - sign(boxlen_d(2),dy)
         if (abs(dz) > boxlen2_d(3)) dz = dz - sign(boxlen_d(3),dz)
         if (abs(dx) + abs(dy) + abs(dz) > ThreeHalf_d*boxlen2_d(1)) then
            dx = dx - sign(boxlen2_d(1),dx)
            dy = dy - sign(boxlen2_d(2),dy)
            dz = dz - sign(boxlen2_d(3),dz)
         end if
      end if
   end if
   r2 = dx**2+dy**2+dz**2

end subroutine PBCr2_cuda

 attributes(device) subroutine PBC_cuda(dx, dy, dz)

 use precision_m
 implicit none
 real(fp_kind), intent(inout)  :: dx, dy, dz

 if (lpbc_d) then
    if (abs(dx) > boxlen2_d(1)) dx = dx - dpbc_d(1) * ANINT(dx*boxleni_d(1))
    if (abs(dy) > boxlen2_d(2)) dy = dy - dpbc_d(2) * ANINT(dy*boxleni_d(2))
    if (abs(dz) > boxlen2_d(3)) dz = dz - dpbc_d(3) * ANINT(dz*boxleni_d(3))
 end if

 end subroutine PBC_cuda


subroutine AllocateDeviceParams


        !use NListModule
        use Random_Module
        implicit none

        integer(4) :: istat

   if(ltime) call CpuAdd('start', 'allocation', 1, uout)
        allocate(iptpt_d(npt,npt))
        !allocate(jpnlist_d(maxnneigh,npartperproc))
        allocate(utwob_d(0:nptpt))
        allocate(ro_d(3,np_alloc))
        allocate(ro_aux(3,np_alloc))
        allocate(r2umin_d(natat))
        allocate(r2atat_d(natat))
        allocate(iubuflow_d(natat))
        !allocate(nneighpn_d(np_alloc))
        write(*,*) "1"
        allocate(iptpn_d(np_alloc))
        allocate(iptpn_aux_d(np_alloc))
        allocate(ubuf_d(nbuf))
        allocate(rotm_d(3,np_alloc))
        allocate(lptm_d(np_alloc))
        allocate(ipnptm_d(np_alloc))
        allocate(dutwob_d(0:nptpt))
        allocate(utwobnew_d(0:nptpt))
        allocate(utwobold_d(0:nptpt))
        allocate(dutwobold(0:nptpt))
        allocate(ucoff_d(natat))
        write(*,*) "2"
        write(*,*) "1"
        allocate(seedsnp(np_alloc))
        write(*,*) "2"
        allocate(seeds_d(np_alloc))
        write(*,*) "3"
        allocate(bondnn_d(2,np))
        allocate(ictpn_d(np_alloc))
        allocate(bondcl_d(4,np_alloc))
        write(*,*) "4"
        allocate(rsumrad(npt,npt))
        write(*,*) "5"
        allocate(rsumrad_h(npt,npt))
        write(*,*) "6"
        allocate(dtran_d(npt))
        write(*,*) "7"
        if (lchain) then
           allocate(bond_d_k(nct))
           bond_d_k = 0.0
           allocate(bond_d_eq(nct))
           bond_d_eq = 0.0
           allocate(bond_d_p(nct))
           bond_d_p = 0
        end if
        allocate(nbondcl_d(np_alloc))
        allocate(ix_d(np_alloc))
        allocate(iy_d(np_alloc))
        allocate(am_d(np_alloc))
        allocate(ipGPU_d(np_alloc))
        allocate(laz_d(na_alloc))
        allocate(laz_aux_d(na_alloc))
        allocate(laztm_d(na_alloc))
        allocate(lspart_d(np_alloc))
        allocate(lchargechange_d(np_alloc))
        allocate(pspart_d(npt))
        allocate(pcharge_d(npt))
        allocate(iananweakcharge_d(na_alloc))
        allocate(iananweakcharge_aux_d(na_alloc))
        allocate(lcounterion_d(np_alloc))
        allocate(ispartmove_d(np_alloc))
        allocate(weight_d(np_alloc))
        allocate(weight_laz(np_alloc))
        allocate(weight_nlaz(np_alloc))
        allocate(weightd_laz(np_alloc))
        allocate(weightd_nlaz(np_alloc))
   if(ltime) call CpuAdd('stop', 'allocation', 1, uout)


end subroutine AllocateDeviceParams

subroutine TransferConstantParams

        use Molmodule
        use Random_Module
        use PotentialModule
        implicit none

        integer(4) :: istat, ipt, jpt, ict, ip, ierra

   if(ltime) call CpuAdd('start', 'transferconstant', 1, uout)
        boxlen2_d = boxlen2
        boxlen_d = boxlen
        boxleni_d = boxleni
        dpbc_d = dpbc
        lPBC_d = lPBC
        lbcbox_d = lbcbox
        lbcrd_d = lbcrd
        lbcto_d = lbcto
        rcut2_d = rcut2
        scrlen_d = scrlen
        nptpt_d = nptpt
        dtran_d = dtran
        r2atat_d = r2atat
        r2umin_d = r2umin
        iubuflow_d = iubuflow
        iptpn_d = iptpn
        iptpn_aux_d = iptpn
        iptpt_d = iptpt
        nbuf_d = nbuf
        ubuf_d = ubuf
        lmonoatom_d = lmonoatom
        lmc_d = lmc
        np_d = np
        npt_d = npt
        lellipsoid_d = lellipsoid
        lsuperball_d = lsuperball
        lptmdutwob_d = lptmdutwob
        iinteractions = np*(np-1)/2
        iinteractions_d = iinteractions
        lchain_d = lchain
        ictpn_d = ictpn
        lcharge_d = lcharge
        lweakcharge_d = lweakcharge
        if (lweakcharge) iananweakcharge_aux_d = iananweakcharge
        if (lweakcharge) lcounterion_d = lcounterion
        pspart_d = pspart
        pcharge_d = pcharge
        if (lchain) then
           do ict =1, nct
              bond_d_k(ict) = bond(ict)%k
              bond_d_eq(ict) = bond(ict)%eq
              bond_d_p(ict) = bond(ict)%p
           end do
           bondnn_d = bondnn
        end if

        lclink_d = lclink
   if(lclink) then
        clink_d_k = clink%k
        clink_d_eq = clink%eq
        clink_d_p = clink%p
        bondcl_d = bondcl
        nbondcl_d = nbondcl
   end if
        lcuda = .true.
        lseq = .true.
        lcuda_mcpass = .false.
        ltest_cuda = .true.

        ro_d = ro
        ro_aux = ro
        if (lweakcharge) laz_d = laz
        if (lweakcharge) laz_aux_d = laz
        sizeofblocks_d = 512
        threadssum =16
        threadssum_d = threadssum

        seeds_d = seedsnp
        beta_d = beta
        do ipt =1, npt
           do jpt = 1, npt
              rsumrad_h(ipt,jpt) = radat(ipt) + radat(jpt)
              rsumrad_h(ipt,jpt) = rsumrad_h(ipt,jpt)**2
           end do
        end do
        do ipt = 1, natat
           ucoff_d(ipt) = ucoff(1,ipt)
        end do
        rsumrad = rsumrad_h
   if(ltime) call CpuAdd('stop', 'transferconstant', 1, uout)
        am_dev = am
        iseed_d = -1228
        iseed_trial_d = iseed_trial
        ipGPU_d = ipGPU
        call TransferCoordinatesToDevice<<<iblock1,256>>>
                  ierra = cudaDeviceSynchronize()
        !call TransferToAux<<<iblock1,256>>>
                  ierra = cudaDeviceSynchronize()
        if(lweakcharge_d) laztm_d = laz_d
        print *, "8"

end subroutine TransferConstantParams

attributes(global) subroutine TransferCoordinatesToDevice

   implicit none
   integer(4) :: id

   id = (blockidx%x-1)*blocksize + threadIDx%x
   if (id <= np_d ) then
      ro_d(1,ipGPU_d(id)) = ro_aux(1,id)
      ro_d(2,ipGPU_d(id)) = ro_aux(2,id)
      ro_d(3,ipGPU_d(id)) = ro_aux(3,id)
      if(lweakcharge_d) laz_d(ipGPU_d(id)) = laz_aux_d(id)
      iptpn_d(ipGPU_d(id)) = iptpn_aux_d(id)
      iananweakcharge_d(ipGPU_d(id)) = ipGPU_d(iananweakcharge_aux_d(id))
   end if


end subroutine TransferCoordinatesToDevice


subroutine TransferDUTotalVarToDevice

        use Molmodule
        implicit none
        integer(4) :: i

        nptm_d = nptm
        ipnptm_d(1:nptm) = ipnptm(1:nptm)
        lptm_d = lptm
        rotm_d(1:3,1:nptm) = rotm(1:3,1:nptm)
        utwobnew_d(0:nptpt) = Zero

end subroutine TransferDUTotalVarToDevice

subroutine CalcWeight_h

   implicit none
   integer(4) :: ipt, xsign

   do ipt = 1, npt
      xsign = sign(One,zat(ipt))
      weight_laz(ipt) = exp(xsign*log(10.0d0)*pHmpK(ipt))
      xsign = -xsign
      weight_nlaz(ipt) = exp(xsign*log(10.0d0)*pHmpK(ipt))
   end do
   weightd_laz = weight_laz
   weightd_nlaz = weight_nlaz

end subroutine


      attributes(device) function Random_dev2(idum)
         implicit none
         integer(k4b), intent(inout) :: idum
         real(8) :: Random_dev2
         integer(k4b), parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
         integer(k4b)   :: k
         !icounter2_d = icounter2_d + 1
         !write(*,*) "icounter2: ", icounter2_d
         if (idum <= 0 .or. iy_dev2 < 0) then           !initialize.
            am_dev=nearest(1.0,-1.0)/im
            iy_dev2=ior(ieor(888889999,abs(idum)),1)
            ix_dev2=ieor(777755555,abs(idum))
            idum=abs(idum)+1                          !set idum positive.
         end if
         ix_dev2=ieor(ix_dev2,ishft(ix_dev2,13))                  !marsaglia shift sequence with period 2^32 − 1.
         ix_dev2=ieor(ix_dev2,ishft(ix_dev2,-17))
         ix_dev2=ieor(ix_dev2,ishft(ix_dev2,5))
         k=iy_dev2/iq                                   !park-miller sequence by schrage’s method, period 2^31 − 2.
         iy_dev2=ia*(iy_dev2-k*iq)-ir*k
         if (iy_dev2 < 0) iy_dev2=iy_dev2+im
         Random_dev2=am_dev*ior(iand(im,ieor(ix_dev2,iy_dev2)),1)     !combine the two generators with masking to ensure nonzero value.
        ! print *, "am ", am_dev
        ! print *, "ix ", ix_dev2
        ! print *, "iy ", iy_dev2
      end function Random_dev2

      attributes(device) function Random_dev(idum)
         implicit none
         integer(k4b), intent(inout) :: idum
         real(8) :: Random_dev
         integer(k4b), parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
         integer(k4b)   :: k
         if (idum <= 0 .or. iy_dev < 0) then           !initialize.
            am_dev=nearest(1.0,-1.0)/im
            iy_dev=ior(ieor(888889999,abs(idum)),1)
            ix_dev=ieor(777755555,abs(idum))
            idum=abs(idum)+1                          !set idum positive.
         end if
         ix_dev=ieor(ix_dev,ishft(ix_dev,13))                  !marsaglia shift sequence with period 2^32 − 1.
         ix_dev=ieor(ix_dev,ishft(ix_dev,-17))
         ix_dev=ieor(ix_dev,ishft(ix_dev,5))
         k=iy_dev/iq                                   !park-miller sequence by schrage’s method, period 2^31 − 2.
         iy_dev=ia*(iy_dev-k*iq)-ir*k
         if (iy_dev < 0) iy_dev=iy_dev+im
         Random_dev=am_dev*ior(iand(im,ieor(ix_dev,iy_dev)),1)     !combine the two generators with masking to ensure nonzero value.
      end function Random_dev


      attributes(device) subroutine Random_d(idum,prandom,id)
         use precision_m
         implicit none
         integer(4), intent(inout) :: idum
         real(fp_kind), intent(inout) :: prandom
         integer(k4b), parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
         integer(4), intent(in) :: id
         integer(k4b) :: k_d
         if (idum <= 0 .or. iy_d(id) < 0) then           !initialize.
            am_d(id)=nearest(1.0,-1.0)/im
            iy_d(id)=ior(ieor(888889999,abs(idum)),1)
            ix_d(id)=ieor(777755555,abs(idum))
            idum=abs(idum)+1                          !set idum positive.
         end if
         ix_d(id)=ieor(ix_d(id),ishft(ix_d(id),13))                  !marsaglia shift sequence with period 2^32 − 1.
         ix_d(id)=ieor(ix_d(id),ishft(ix_d(id),-17))
         ix_d(id)=ieor(ix_d(id),ishft(ix_d(id),5))
         k_d=iy_d(id)/iq                                   !park-miller sequence by schrage’s method, period 2^31 − 2.
         iy_d(id)=ia*(iy_d(id)-k_d*iq)-ir*k_d
         if (iy_d(id) < 0) iy_d(id)=iy_d(id)+im
         prandom= am_d(id)*ior(iand(im,ieor(ix_d(id),iy_d(id))),1)     !combine the two generators with masking to ensure nonzero value.
      end subroutine Random_d

      subroutine GenerateSeeds

         use Random_Module
         implicit none
         integer(4) :: ip
         ix_dev = ix
         iy_dev = iy
         am_dev = am
         iseed_d = -1228
         do ip = 1, np
           seedsnp(ip) = Random_int(iseed)
           seedsnp(ip) = -abs(seedsnp(ip))
         end do
         seeds_d = seedsnp

      end subroutine

      function Random_int(idum)
         use Random_Module
         implicit none
         integer(k4b), intent(inout) :: idum
         integer(k4b) :: Random_int
         integer(k4b), parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
         integer(k4b)   :: k
         integer(k4b)     :: a = 1000, b = 100000
         if (idum <= 0 .or. iy < 0) then           !initialize.
            am=nearest(1.0,-1.0)/im
            iy=ior(ieor(888889999,abs(idum)),1)
            ix=ieor(777755555,abs(idum))
            idum=abs(idum)+1                          !set idum positive.
         end if
         ix=ieor(ix,ishft(ix,13))                  !marsaglia shift sequence with period 2^32 − 1.
         ix=ieor(ix,ishft(ix,-17))
         ix=ieor(ix,ishft(ix,5))
         k=iy/iq                                   !park-miller sequence by schrage’s method, period 2^31 − 2.
         iy=ia*(iy-k*iq)-ir*k
         if (iy < 0) iy=iy+im
         Random_int = ceiling(a+(b-a)*am*ior(iand(im,ieor(ix,iy)),1))     !combine the two generators with masking to ensure nonzero value.
      end function Random_int

      attributes(global) subroutine ConvertCoordinatesAuxToDevice

         implicit none
         integer(4) :: id

         id = (blockidx%x-1)*blocksize + threadIDx%x
         if (id <= np_d) then
            ro_d(1,id) = ro_aux(1,id)
            ro_d(2,id) = ro_aux(2,id)
            ro_d(3,id) = ro_aux(3,id)
            if (lweakcharge_d) laz_d(id) = laz_aux_d(id)
         end if

      end subroutine ConvertCoordinatesAuxToDevice

      attributes(global) subroutine ConvertCoordinates

         implicit none
         integer(4) :: id

         id = (blockidx%x-1)*blocksize + threadIDx%x
         if (id <= np_d) then
            ro_d(1,id) = ro_aux(1,ipGPU_d(id))
            ro_d(2,id) = ro_aux(2,ipGPU_d(id))
            ro_d(3,id) = ro_aux(3,ipGPU_d(id))
            if (lweakcharge_d) laz_d(id)  = laz_aux_d(ipGPU_d(id))
         end if

      end subroutine ConvertCoordinates


end module CUDAModule
