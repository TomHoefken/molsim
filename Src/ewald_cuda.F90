module EwaldCudaModule

      use CUDAModule
      use EnergyModule
      implicit none


      contains

attributes(global) subroutine DUTwoBodyEwaldRecStd_cuda

   use EnergyModule
   #if defined (_CUDA_)
      use CUDAModule
   #endif
   implicit none
   !character(40), parameter :: txroutine ='UEwaldRecStd'
   integer(4) :: kn, nx, ny, nz, ia, ialoc, ikvec2
   real(8)    :: term, termnew, termold

!   if (ltime) call CpuAdd('start', txroutine, 3, uout)

! ... calculate eikxtm, eikytm, and eikztm for moving particles

   call EwaldSetArrayTM_cuda

   kn = kvecoffmyid_d
   ikvec2 = 0
   do nz = 0, ncut_d
      do ny = 0, ncut_d
         if (ny**2+nz**2 > ncut2_d) cycle
         ikvec2 = ikvec2+1
         !print *, ikvec2, kvecmyid_d(1), kvecmyid_d(2)
         if (ikvec2 < kvecmyid_d(1) .or. ikvec2 > kvecmyid_d(2)) cycle  ! parallelize over k-vectors
         do ialoc = 1, natm_d
            ia = ianatm_d(ialoc)
            eikyzm_d(ia)      = conjg(eiky_d(ia,ny))     *eikz_d(ia,nz)
            eikyzp_d(ia)      =       eiky_d(ia,ny)      *eikz_d(ia,nz)
            eikyzmtm_d(ialoc) = conjg(eikytm_d(ialoc,ny))*eikztm_d(ialoc,nz)
            eikyzptm_d(ialoc) =       eikytm_d(ialoc,ny) *eikztm_d(ialoc,nz)
         end do
         !print *, "eikyzm_d: ", eikyzm_d(1)

         do nx = 0, ncut_d
            if ((lbcrd_d .or. lbcto_d) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2_d) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn + 1
            sumeikrtm_d(kn,1) = sumeikr_d(kn,1)
            sumeikrtm_d(kn,2) = sumeikr_d(kn,2)
            sumeikrtm_d(kn,3) = sumeikr_d(kn,3)
            sumeikrtm_d(kn,4) = sumeikr_d(kn,4)
            do ialoc = 1, natm_d
               ia = ianatm_d(ialoc)
               sumeikrtm_d(kn,1) = sumeikrtm_d(kn,1)+az_d(ia)*  &
                  (conjg(eikxtm_d(ialoc,nx))*eikyzmtm_d(ialoc) - conjg(eikx_d(ia,nx))*eikyzm_d(ia))
               sumeikrtm_d(kn,2) = sumeikrtm_d(kn,2)+az_d(ia)*  &
                  (conjg(eikxtm_d(ialoc,nx))*eikyzptm_d(ialoc) - conjg(eikx_d(ia,nx))*eikyzp_d(ia))
               sumeikrtm_d(kn,3) = sumeikrtm_d(kn,3)+az_d(ia)*  &
                        (eikxtm_d(ialoc,nx) *eikyzmtm_d(ialoc) -       eikx_d(ia,nx) *eikyzm_d(ia))
               sumeikrtm_d(kn,4) = sumeikrtm_d(kn,4)+az_d(ia)*  &
                        (eikxtm_d(ialoc,nx) *eikyzptm_d(ialoc) -       eikx_d(ia,nx) *eikyzp_d(ia))
            end do

            termnew = real(sumeikrtm_d(kn,1))**2 + aimag(sumeikrtm_d(kn,1))**2 + real(sumeikrtm_d(kn,2))**2 + aimag(sumeikrtm_d(kn,2))**2 &
                    + real(sumeikrtm_d(kn,3))**2 + aimag(sumeikrtm_d(kn,3))**2 + real(sumeikrtm_d(kn,4))**2 + aimag(sumeikrtm_d(kn,4))**2
            termold = real(sumeikr_d(kn,1))**2   + aimag(sumeikr_d(kn,1))**2   + real(sumeikr_d(kn,2))**2   + aimag(sumeikr_d(kn,2))**2 &
                    + real(sumeikr_d(kn,3))**2   + aimag(sumeikr_d(kn,3))**2   + real(sumeikr_d(kn,4))**2   + aimag(sumeikr_d(kn,4))**2
            term    = kfac_d(kn)*(termnew - termold)
            durec_d   = durec_d + term

         end do
      end do
   end do
   !print *, kvecmyid(1), kvecmyid(2), kvecoffmyid
   !print *, ncut, ikvec2, kn

!   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

end subroutine DUTwoBodyEwaldRecStd_cuda

attributes(global) subroutine DUTwoBodyEwaldRecStd_cuda_par

   use EnergyModule
   #if defined (_CUDA_)
      use CUDAModule
   #endif
   implicit none
   !character(40), parameter :: txroutine ='UEwaldRecStd'
   integer(4) :: kn, nx, ny, nz!, ia, ialoc, id_int
   integer(4),shared :: ia, ialoc, id_int
   real(8),shared :: term(256), termnew(256), termold(256)

!   if (ltime) call CpuAdd('start', txroutine, 3, uout)
   kn = ((blockIDx%x-1) * blocksize + threadIDx%x)
   nx = kfacnx_d(kn)
   ny = kfacny_d(kn)
   nz = kfacnz_d(kn)
   id_int = threadIDx%x
   term(kn) = 0.0

! ... calculate eikxtm, eikytm, and eikztm for moving particles

   if (kn == 1) call EwaldSetArrayTM_cuda
   call syncthreads

      if (kn <= nkvec_d) then
         do ialoc = 1, natm_d
            ia = ianatm_d(ialoc)
            eikyzm2_d(kn,ia)      = conjg(eiky_d(ia,ny))     *eikz_d(ia,nz)
            eikyzp2_d(kn,ia)      =       eiky_d(ia,ny)      *eikz_d(ia,nz)
            eikyzmtm2_d(kn,ialoc) = conjg(eikytm_d(ialoc,ny))*eikztm_d(ialoc,nz)
            eikyzptm2_d(kn,ialoc) =       eikytm_d(ialoc,ny) *eikztm_d(ialoc,nz)
         end do

            sumeikrtm_d(kn,1) = sumeikr_d(kn,1)
            sumeikrtm_d(kn,2) = sumeikr_d(kn,2)
            sumeikrtm_d(kn,3) = sumeikr_d(kn,3)
            sumeikrtm_d(kn,4) = sumeikr_d(kn,4)
            do ialoc = 1, natm_d
               ia = ianatm_d(ialoc)
               sumeikrtm_d(kn,1) = sumeikrtm_d(kn,1)+az_d(ia)*  &
                  (conjg(eikxtm_d(ialoc,nx))*eikyzmtm2_d(kn,ialoc) - conjg(eikx_d(ia,nx))*eikyzm2_d(kn,ia))
               sumeikrtm_d(kn,2) = sumeikrtm_d(kn,2)+az_d(ia)*  &
                  (conjg(eikxtm_d(ialoc,nx))*eikyzptm2_d(kn,ialoc) - conjg(eikx_d(ia,nx))*eikyzp2_d(kn,ia))
               sumeikrtm_d(kn,3) = sumeikrtm_d(kn,3)+az_d(ia)*  &
                        (eikxtm_d(ialoc,nx) *eikyzmtm2_d(kn,ialoc) -       eikx_d(ia,nx) *eikyzm2_d(kn,ia))
               sumeikrtm_d(kn,4) = sumeikrtm_d(kn,4)+az_d(ia)*  &
                        (eikxtm_d(ialoc,nx) *eikyzptm2_d(kn,ialoc) -       eikx_d(ia,nx) *eikyzp2_d(kn,ia))
            end do

            termnew(kn) = real(sumeikrtm_d(kn,1))**2 + aimag(sumeikrtm_d(kn,1))**2 + real(sumeikrtm_d(kn,2))**2 + aimag(sumeikrtm_d(kn,2))**2 &
                    + real(sumeikrtm_d(kn,3))**2 + aimag(sumeikrtm_d(kn,3))**2 + real(sumeikrtm_d(kn,4))**2 + aimag(sumeikrtm_d(kn,4))**2
            termold(kn) = real(sumeikr_d(kn,1))**2   + aimag(sumeikr_d(kn,1))**2   + real(sumeikr_d(kn,2))**2   + aimag(sumeikr_d(kn,2))**2 &
                    + real(sumeikr_d(kn,3))**2   + aimag(sumeikr_d(kn,3))**2   + real(sumeikr_d(kn,4))**2   + aimag(sumeikr_d(kn,4))**2
            term(kn)    = kfac_d(kn)*(termnew(kn) - termold(kn))
      end if
      termold(kn) = AtomicAdd(durec_d,term(kn))

!   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

end subroutine DUTwoBodyEwaldRecStd_cuda_par

attributes(device) subroutine EwaldSetArrayTM_cuda

   !use EnergyModule
   implicit none

   integer(4) :: ialoc, icut

   do ialoc = 1, natm_d
      eikxtm_d(ialoc,0) = cmplx(One_d,Zero_d)
      eikytm_d(ialoc,0) = cmplx(One_d,Zero_d)
      eikztm_d(ialoc,0) = cmplx(One_d,Zero_d)
      eikxtm_d(ialoc,1) = cmplx(cos(TwoPiBoxi_d(1)*rtm_d(1,ialoc)),sin(TwoPiBoxi_d(1)*rtm_d(1,ialoc)))
      eikytm_d(ialoc,1) = cmplx(cos(TwoPiBoxi_d(2)*rtm_d(2,ialoc)),sin(TwoPiBoxi_d(2)*rtm_d(2,ialoc)))
      eikztm_d(ialoc,1) = cmplx(cos(TwoPiBoxi_d(3)*rtm_d(3,ialoc)),sin(TwoPiBoxi_d(3)*rtm_d(3,ialoc)))
   end do
   do icut = 2, ncut_d
      do ialoc = 1, natm_d
         eikxtm_d(ialoc,icut) = eikxtm_d(ialoc,icut-1)*eikxtm_d(ialoc,1)
         eikytm_d(ialoc,icut) = eikytm_d(ialoc,icut-1)*eikytm_d(ialoc,1)
         eikztm_d(ialoc,icut) = eikztm_d(ialoc,icut-1)*eikztm_d(ialoc,1)
      end do
   end do

end subroutine EwaldSetArrayTM_cuda



end module EwaldCudaModule
