module EwaldCudaModule

      use mol_cuda
      !use EnergyModule
      implicit none


   ! for Ewald Summation
      contains

attributes(global) subroutine DUTwoBodyEwaldRecStd_cuda

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

   implicit none
   !character(40), parameter :: txroutine ='UEwaldRecStd'
   integer(4) :: kn, nx, ny, nz!, ia, ialoc, id_int
   integer(4),shared :: ia, ialoc, id_int
   real(8),shared :: term(256), termnew(256), termold(256)
   integer(4)   :: blocksize = 256

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
   real(8) :: One_d = 1.0d0
   real(8) :: Zero_d = 0.0d0

   do ialoc = 1, natm_d
      eikxtm_d(ialoc,0) = cmplx(One_d,Zero_d)
      eikytm_d(ialoc,0) = cmplx(One_d,Zero_d)
      eikztm_d(ialoc,0) = cmplx(One_d,Zero_d)
      eikxtm_d(ialoc,1) = cmplx(cos(TwoPiBoxi_d(1)*rotm_d(1,ialoc)),sin(TwoPiBoxi_d(1)*rotm_d(1,ialoc)))
      eikytm_d(ialoc,1) = cmplx(cos(TwoPiBoxi_d(2)*rotm_d(2,ialoc)),sin(TwoPiBoxi_d(2)*rotm_d(2,ialoc)))
      eikztm_d(ialoc,1) = cmplx(cos(TwoPiBoxi_d(3)*rotm_d(3,ialoc)),sin(TwoPiBoxi_d(3)*rotm_d(3,ialoc)))
   end do
   do icut = 2, ncut_d
      do ialoc = 1, natm_d
         eikxtm_d(ialoc,icut) = eikxtm_d(ialoc,icut-1)*eikxtm_d(ialoc,1)
         eikytm_d(ialoc,icut) = eikytm_d(ialoc,icut-1)*eikytm_d(ialoc,1)
         eikztm_d(ialoc,icut) = eikztm_d(ialoc,icut-1)*eikztm_d(ialoc,1)
      end do
   end do

end subroutine EwaldSetArrayTM_cuda

attributes(device) subroutine EwaldUpdateArray_cuda(ia)

   implicit none
   integer(4), intent(in) :: ia
   integer(4) :: icut,ialoc

      do icut = 0, ncut_d
         do ialoc = 1, natm_d
            eikx_d(ia,icut) = eikxtm_d(ialoc,icut)
            eiky_d(ia,icut) = eikytm_d(ialoc,icut)
            eikz_d(ia,icut) = eikztm_d(ialoc,icut)
         end do
      end do
      sumeikr_d(1:nkvec_d,1:4) = sumeikrtm_d(1:nkvec_d,1:4)
end subroutine EwaldUpdateArray_cuda

subroutine AllocateEwaldAuxParams

   use EnergyModule
   implicit none
   integer(4) :: istat, imem

   !imem = ncut+1
   !allocate(eikx_aux(na,imem))
   !allocate(eiky_aux(na,imem))
   !allocate(eikz_aux(na,imem))
   !istat = cudaMallocPitch(eikx_aux(na,0:ncut),na,na,imem)
   !istat = cudaMallocPitch(eiky_aux(na,0:ncut),na,na,imem)
   !istat = cudaMallocPitch(eikz_aux(na,0:ncut),na,na,imem)

end subroutine AllocateEwaldAuxParams

subroutine TransferEwaldAuxParams

   use EnergyModule
   implicit none

   eikx_aux = eikx
   eiky_aux = eiky
   eikz_aux = eikz
   eikx_d = eikx
   eiky_d = eiky
   eikz_d = eikz

end subroutine TransferEwaldAuxParams

attributes(device) subroutine TransferEwaldtoDevice(id)

   !use EwaldCUDAModule
   implicit none
   integer(4), intent(in)  :: id
   integer(4) :: icut

   do icut = 0, ncut_d
      eikx_d(ipGPU_d(id),icut) = eikx_aux(id,icut)
      eiky_d(ipGPU_d(id),icut) = eiky_aux(id,icut)
      eikz_d(ipGPU_d(id),icut) = eikz_aux(id,icut)
   end do


end subroutine TransferEwaldtoDevice

attributes(device) subroutine DUTwoBodyEwaldRecStd_liang(ip)

   implicit none
   !character(40), parameter :: txroutine ='UEwaldRecStd'
   integer(4) :: kn, nx, ny, nz!, ia, ialoc, id_int
   integer(4) :: ia, ialoc, id_int
   real(8)    :: term(256), termnew(256), termold(256)
   integer(4)   :: blocksize = 256
   integer(4), intent(inout) :: ip

!   if (ltime) call CpuAdd('start', txroutine, 3, uout)
   kn = ((blockIDx%x-1) * blocksize + threadIDx%x)
   id_int = threadIDx%x
   natm_d = 1

! ... calculate eikxtm, eikytm, and eikztm for moving particles

   !if (kn == 1) call EwaldSetArrayTM_liang(ip)
   call syncthreads

      if (kn <= nkvec_d) then
             term(kn) = 0.0
            nx = kfacnx_d(kn)
            ny = kfacny_d(kn)
            nz = kfacnz_d(kn)
         !do ialoc = 1, natm_d
            !ia = ianatm_d(ialoc)
            ia = ip
            eikyzm2_d(kn,ia)      = conjg(eiky_d(ia,ny))     *eikz_d(ia,nz)
            eikyzp2_d(kn,ia)      =       eiky_d(ia,ny)      *eikz_d(ia,nz)
            eikyzmtm2_d(kn,ialoc) = conjg(eikytm_d(ialoc,ny))*eikztm_d(ialoc,nz)
            eikyzptm2_d(kn,ialoc) =       eikytm_d(ialoc,ny) *eikztm_d(ialoc,nz)
         !end do

            sumeikrtm_d(kn,1) = sumeikr_d(kn,1)
            sumeikrtm_d(kn,2) = sumeikr_d(kn,2)
            sumeikrtm_d(kn,3) = sumeikr_d(kn,3)
            sumeikrtm_d(kn,4) = sumeikr_d(kn,4)
            !do ialoc = 1, natm_d
               !ia = ianatm_d(ialoc)
               ia = ip
               sumeikrtm_d(kn,1) = sumeikrtm_d(kn,1)+az_d(ia)*  &
                  (conjg(eikxtm_d(ialoc,nx))*eikyzmtm2_d(kn,ialoc) - conjg(eikx_d(ia,nx))*eikyzm2_d(kn,ia))
               sumeikrtm_d(kn,2) = sumeikrtm_d(kn,2)+az_d(ia)*  &
                  (conjg(eikxtm_d(ialoc,nx))*eikyzptm2_d(kn,ialoc) - conjg(eikx_d(ia,nx))*eikyzp2_d(kn,ia))
               sumeikrtm_d(kn,3) = sumeikrtm_d(kn,3)+az_d(ia)*  &
                        (eikxtm_d(ialoc,nx) *eikyzmtm2_d(kn,ialoc) -       eikx_d(ia,nx) *eikyzm2_d(kn,ia))
               sumeikrtm_d(kn,4) = sumeikrtm_d(kn,4)+az_d(ia)*  &
                        (eikxtm_d(ialoc,nx) *eikyzptm2_d(kn,ialoc) -       eikx_d(ia,nx) *eikyzp2_d(kn,ia))
            !end do

            termnew(kn) = real(sumeikrtm_d(kn,1))**2 + aimag(sumeikrtm_d(kn,1))**2 + real(sumeikrtm_d(kn,2))**2 + aimag(sumeikrtm_d(kn,2))**2 &
                    + real(sumeikrtm_d(kn,3))**2 + aimag(sumeikrtm_d(kn,3))**2 + real(sumeikrtm_d(kn,4))**2 + aimag(sumeikrtm_d(kn,4))**2
            termold(kn) = real(sumeikr_d(kn,1))**2   + aimag(sumeikr_d(kn,1))**2   + real(sumeikr_d(kn,2))**2   + aimag(sumeikr_d(kn,2))**2 &
                    + real(sumeikr_d(kn,3))**2   + aimag(sumeikr_d(kn,3))**2   + real(sumeikr_d(kn,4))**2   + aimag(sumeikr_d(kn,4))**2
            term(kn)    = kfac_d(kn)*(termnew(kn) - termold(kn))
            !termold(kn) = AtomicAdd(durec_d,term(kn))
      end if

!   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

end subroutine DUTwoBodyEwaldRecStd_liang

attributes(device) subroutine DUTwoBodyEwaldRecStd_liang_ser(ip)

   implicit none
   !character(40), parameter :: txroutine ='UEwaldRecStd'
   integer(4) :: kn, nx, ny, nz, ia, ialoc, ikvec2
   real(8)    :: term, termnew, termold
   integer(4), intent(inout) :: ip

!   if (ltime) call CpuAdd('start', txroutine, 3, uout)

! ... calculate eikxtm, eikytm, and eikztm for moving particles

   natm_d = 1
   call EwaldSetArrayTM_liang(ip)

   kn = kvecoffmyid_d
   ikvec2 = 0
   do nz = 0, ncut_d
      do ny = 0, ncut_d
         if (ny**2+nz**2 > ncut2_d) cycle
         ikvec2 = ikvec2+1
         !print *, ikvec2, kvecmyid_d(1), kvecmyid_d(2)
         if (ikvec2 < kvecmyid_d(1) .or. ikvec2 > kvecmyid_d(2)) cycle  ! parallelize over k-vectors
         do ialoc = 1, natm_d
            ia = ip
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
               ia = ip
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

end subroutine DUTwoBodyEwaldRecStd_liang_ser

attributes(device) subroutine EwaldSetArrayTM_liang(ip)

   !use EnergyModule
   implicit none

   integer(4) :: ialoc, icut
   real(8) :: One_d = 1.0d0
   real(8) :: Zero_d = 0.0d0
   integer(4), intent(in) :: ip

   natm_d = 1
   do ialoc = 1, natm_d
      eikxtm_d(ialoc,0) = cmplx(One_d,Zero_d)
      eikytm_d(ialoc,0) = cmplx(One_d,Zero_d)
      eikztm_d(ialoc,0) = cmplx(One_d,Zero_d)
      eikxtm_d(ialoc,1) = cmplx(cos(TwoPiBoxi_d(1)*rotm_d(1,ip)),sin(TwoPiBoxi_d(1)*rotm_d(1,ip)))
      eikytm_d(ialoc,1) = cmplx(cos(TwoPiBoxi_d(2)*rotm_d(2,ip)),sin(TwoPiBoxi_d(2)*rotm_d(2,ip)))
      eikztm_d(ialoc,1) = cmplx(cos(TwoPiBoxi_d(3)*rotm_d(3,ip)),sin(TwoPiBoxi_d(3)*rotm_d(3,ip)))
   end do
   do icut = 2, ncut_d
      do ialoc = 1, natm_d
         eikxtm_d(ialoc,icut) = eikxtm_d(ialoc,icut-1)*eikxtm_d(ialoc,1)
         eikytm_d(ialoc,icut) = eikytm_d(ialoc,icut-1)*eikytm_d(ialoc,1)
         eikztm_d(ialoc,icut) = eikztm_d(ialoc,icut-1)*eikztm_d(ialoc,1)
      end do
   end do

end subroutine EwaldSetArrayTM_liang

attributes(device) subroutine EwaldUpdateArray_liang(ia)

   implicit none
   integer(4), intent(in) :: ia
   integer(4) :: icut,ialoc

   natm_d = 1
      do icut = 0, ncut_d
         do ialoc = 1, natm_d
            eikx_d(ia,icut) = eikxtm_d(ialoc,icut)
            eiky_d(ia,icut) = eikytm_d(ialoc,icut)
            eikz_d(ia,icut) = eikztm_d(ialoc,icut)
         end do
      end do
      sumeikr_d(1:nkvec_d,1:4) = sumeikrtm_d(1:nkvec_d,1:4)

end subroutine EwaldUpdateArray_liang
end module EwaldCudaModule
