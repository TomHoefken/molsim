module mol_cuda

    use Molmodule
    use cudafor
    use precision_m
    implicit none
   logical :: lcuda = .true.

   real(fp_kind),    constant :: ThreeHalf_d = 1.5
   real(fp_kind),    constant :: SqTwo_d       = sqrt(Two)
   logical,device       :: lbcbox_d                 ! box-like cell (rätblock)
   logical,device       :: lbcrd_d                  ! rhombic dodecahedral cell
   logical,device       :: lbcto_d                  ! truncated octahedral cell
   real(fp_kind),device       :: boxlen_d(3)
   real(fp_kind),device       :: boxlen2_d(3)             ! boxlen/2
   real(fp_kind),device       :: boxleni_d(3)             ! boxlen/2
   real(fp_kind),device       :: dpbc_d(3)             ! /2
   logical,device       :: lPBC_d                   ! periodic boundary conditions
   integer(4),device              :: np_d           ! number of particles
   integer(4), device    :: nptpt_d                  ! number of different particle type pairs
   integer(4), device    :: npt_d
   logical,device                    :: lmonoatom_d
   real(fp_kind), device, allocatable       :: r2atat_d(:)     !
   integer(4), device, allocatable :: iptpn_d(:)     ! particle (1:np)               -> its particle type (1:npt)
   integer(4), device, allocatable :: iptpn_aux_d(:)     ! particle (1:np)               -> its particle type (1:npt)

   integer(4),device, allocatable :: iptpt_d(:,:)   ! two particle types (1:npt)    -> particle type pair (1:nptpt)

   logical,device       :: lmc_d                    ! flag for monte carlo simulation
   real(fp_kind),device       :: virial_d                 ! virial
   real(fp_kind), device, allocatable :: ro_d(:,:)         ! particle position
   real(fp_kind), device, allocatable :: ro_aux(:,:)         ! particle position
   integer(4),device, allocatable :: nneighpn_d(:)  ! particle (local) -> number of neighbours
   integer(4),device, allocatable :: jpnlist_d(:,:) ! ineigh (local list) and ip (global or local) -> neigbour particle (1:np)
   integer(4),device              :: nbuf_d         ! length of buffer
   real(fp_kind), device,    allocatable :: ubuf_d(:)      ! buffer for potential table
   real(fp_kind), device                 :: rcut2_d        ! rcut**2
   real(fp_kind), device,    allocatable :: r2umin_d(:)    ! lower limit squared of potential table
   integer(4),device, allocatable :: iubuflow_d(:)  ! points on the first entry for iatjat
   real(fp_kind), device,allocatable :: utwob_d(:)
   real(fp_kind), device  :: utot_d
   real(fp_kind), device  :: virtwob_d

!... potential
   real(fp_kind), device, allocatable :: ucoff_d(:)
   real(fp_kind), device :: scrlen_d

! MCPass
   integer(4),device           :: iptmove_d
   integer(4),device           :: ipmove_d
   real(fp_kind), device, allocatable :: dtran_d(:)
!... in DuTotal
   logical, device      :: lhsoverlap_d
   integer(4),device    :: nptm_d
   integer(4), device, allocatable :: ipnptm_d(:)
   logical, device, allocatable    :: lptm_d(:)
   real(fp_kind), device, allocatable    :: rotm_d(:,:)
   logical, device                 :: lellipsoid_d
   logical, device                 :: lsuperball_d
   real(fp_kind), device, allocatable    :: dutwob_d(:)
   logical,device       :: lptmdutwob_d             ! flag for calulating dutobdy among moving particles
   real(fp_kind), device,allocatable :: utwobnew_d(:)
   real(fp_kind), device,allocatable :: utwobold_d(:)
   real(fp_kind), allocatable :: dutwobold(:)


   integer(4) :: iinteractions
   integer(4),device :: iinteractions_d
   integer(4), device :: ierror_d
   integer(4), device :: sizeofblocks_d

   integer(4) :: threadssum
   integer(4),device :: threadssum_d

   ! bonds
   real(fp_kind), device, allocatable :: bond_d_k(:)
   real(fp_kind), device, allocatable :: bond_d_eq(:)
   real(fp_kind), device, allocatable :: bond_d_p(:)
   integer(4), device, allocatable :: bondnn_d(:,:)
   real(fp_kind) :: bond_aux
   logical, device :: lchain_d
   integer(4), device, allocatable :: ictpn_d(:)
   real(fp_kind), device, allocatable :: rsumrad(:,:)
   real(fp_kind), device :: clink_d_k
   real(fp_kind), device :: clink_d_eq
   real(fp_kind), device :: clink_d_p
   integer(4), device, allocatable :: bondcl_d(:,:)
   integer(4), device, allocatable :: nbondcl_d(:)
   logical, device :: lclink_d
   real(fp_kind), allocatable :: rsumrad_h(:,:)
   real(fp_kind), device, allocatable :: sig(:)
   real(fp_kind), device, allocatable :: eps(:)

   integer(4), allocatable :: seedsnp(:)

   real(fp_kind), device :: beta_d

   logical :: lseq
   logical :: lcuda_mcpass
   logical,device :: ltest_cuda


   real(8) :: u_aux

      !for random generators
   real(8),device :: am_dev
   integer(k4b),device :: ix_dev=-1,iy_dev=-1
   integer(k4b),device :: ix_dev2=-1,iy_dev2=-1
   integer, constant :: k4b_d=selected_int_kind(9) ! = 4 on intel fortran and gfortran
   real(8), device, allocatable :: am_d(:)
   integer(k4b),device, allocatable :: ix_d(:),iy_d(:)
   integer(4), allocatable :: seeds(:)
   integer(4),device, allocatable :: seeds_d(:)
   integer(4) :: icounter
   integer(4),device :: icounter_d
   integer(4) :: icounter2
   integer(4),device :: icounter2_d
   integer(k4b), device :: iseed_d
   integer(k4b), device :: iseed2_d = -1348
   integer(k4b), device :: iseed_trial_d

   logical, device :: lcharge_d
   logical, device :: lweakcharge_d
   logical, device, allocatable :: laz_d(:)
   logical, device, allocatable :: laz_aux_d(:)
   logical, device, allocatable :: laztm_d(:)
   logical, device, allocatable :: lspart_d(:)
   logical, device, allocatable :: lchargechange_d(:)

   real(fp_kind), device, allocatable :: pspart_d(:)
   real(fp_kind), device, allocatable :: pcharge_d(:)
   integer(4), device, allocatable :: iananweakcharge_d(:)
   integer(4), device, allocatable :: iananweakcharge_aux_d(:)
   logical, device, allocatable :: lcounterion_d(:)


   real(fp_kind), device, allocatable :: weight_d(:)
   real(fp_kind),         allocatable :: weight_laz(:)
   real(fp_kind),         allocatable :: weight_nlaz(:)
   real(fp_kind), device, allocatable :: weightd_laz(:)
   real(fp_kind), device, allocatable :: weightd_nlaz(:)



end module mol_cuda
