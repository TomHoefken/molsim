module NetworkSynModule

   implicit none
   logical, allocatable       :: lactivestart(:)
   logical, allocatable       :: lactive(:)
   real(8), allocatable       :: pbond(:,:)
   logical                    :: lsyn
   real(8), allocatable       :: pactivestart(:)
   integer(4)                 :: iptcl

   contains

      subroutine NetworkSynDriver(iStage)

         use MolModule

         implicit none
         integer(4), intent(in) :: iStage



         select case (iStage)
         case (iReadInput)
            call SetParams

         case (iBeforeSimulation)

            call SetActiveParticles

         case (iSimulationStep)

            call MakeSynthesis3

         end select

      end subroutine NetworkSynDriver

      subroutine SetParams

         use MolModule
         implicit none

         namelist /nmlNetworkSyn/ lactivestart, pbond, pactivestart, iptcl

         allocate(lactivestart(npt))
         allocate(lactive(np))
         allocate(pbond(npt,npt))
         allocate(pactivestart(npt))
         if (.not. allocated(nbondcl)) then
            allocate(nbondcl(np_alloc))
            nbondcl = 0
         end if
         if (.not. allocated(bondcl)) then
            allocate(bondcl(4,np_alloc))
            bondcl = 0
         end if
         if (.not. allocated(bondnn)) then
            allocate(bondnn(2,np_alloc))
            bondnn = 0
         end if
         if (.not. allocated(ictpn)) then
            allocate(ictpn(np_alloc))
            ictpn = 0
         end if
         if (.not. allocated(icnpn)) then
            allocate(icnpn(np_alloc))
            icnpn = 0
         end if
         if (allocated(ictcn)) deallocate(ictcn)
         allocate(ictcn(100))

         lactivestart = .false.
         lactive = .false.
         pbond = 0.0
         pactivestart = 0.0
         iptcl = 0
         npct = 0
         nc = 0
         ictcn = 0
         nct = 0
         ictpn = 0
         icnpn = 0
         bondnn = 0
         bondcl = 0
         nbondcl = 0

         rewind(uin)
         read(uin,nmlNetworkSyn)

      end subroutine SetParams

      subroutine SetActiveParticles

         use Random_Module
         use MolModule
         implicit none
         integer(4) :: ip
         real(8) :: Random2

         do ip = 1, np

            if (.not. lactivestart(iptpn(ip))) cycle
            if (Random2(iseed) < pactivestart(iptpn(ip))) then
               lactive(ip) = .true.
               print *, ip, lactive(ip)
               if (iptpn(ip) /= iptcl) then
                  nc = nc +1
                  ipnsegcn(1,nc) = ip
                  nct = nct + 1
                  ictcn(nc) = nc
                  !print *, "nc: ", nc, ictcn(nc)
                  ictpn(ip) = nc
                  icnpn(ip) = nc
                  npct(nc) = 1
               end if
            end if

         end do
         lchain = .true.
         lclink = .true.
         print *, "number of chain types", nct
         if (.not. allocated(angle)) then
            !allocate(angle(nct))
            allocate(angle(100))
            angle = bond_var(0.0, 0, 0.0)
         end if
         if (.not. allocated(bond)) then
            !allocate(bond(nct))
            allocate(bond(100))
            bond = bond_var(0.1, 2, 5.0)
         end if
            clink = bond_var(1.0, 2, 6.0)
         if (.not. allocated(ipnsegcn)) then
            allocate(ipnsegcn(500,100))
            ipnsegcn = 0
         end if
         !print *, "bondeq: ", bond%eq

      end subroutine SetActiveParticles

      subroutine MakeSynthesis

         use MolModule
         use Random_Module
         implicit none
         integer(4) :: ip, jp, ipt, iseg, qp, icttemp, npcttemp, i
         real(8) :: dx, dy, dz, rdist, Random2
         logical :: llinked

         do ip = 1, np
            if (.not. lactive(ip)) cycle
            if (bondnn(2,ip) /= 0) cycle
            if (nbondcl(ip) == maxnbondcl(iptpn(ip))) cycle
            do jp = 1, np
               !llinked = .false.
               if (.not. lactive(ip)) exit
               if (ip == jp) cycle
               if (bondnn(2,jp) /= 0) cycle
               if (icnpn(ip) == icnpn(jp)) cycle
               !do i = 1, 4
               !   if (bondcl(i,ip) == jp) llinked =.true.
               !end do
               !if (llinked) cycle
               !if (any(bondcl(:,ip) == jp)) cycle
               dx = ro(1,ip) - ro(1,jp)
               dy = ro(2,ip) - ro(2,jp)
               dz = ro(3,ip) - ro(3,jp)
               call PBCr2(dx, dy, dz, rdist)
               if (rdist > 50.0) cycle   !(1.2*bondeq**2)
               print *, ip, jp, rdist
               print *, lactive(ip), lactive(jp)
               if (pbond(iptpn(ip), iptpn(jp)) < Random2(iseed)) cycle
                  ! if one particle is crosslinker
               if (iptcl == iptpn(ip) .or. iptcl == iptpn(jp)) then
                  !if (nbondcl(ip) == maxnbondcl(iptpn(ip))) cycle
                  if (nbondcl(jp) == maxnbondcl(iptpn(jp))) cycle
                 nbondcl(ip) = nbondcl(ip) + 1
                 nbondcl(jp) = nbondcl(jp) + 1
                 bondcl(nbondcl(ip),ip) = jp
                 bondcl(nbondcl(jp),jp) = ip
                 if (iptcl == iptpn(ip)) then
                    if (ictpn(jp) == 0) then
                       nc = nc + 1
                       nct = nct + 1
                       ictcn(nc) = nc
                       ipnsegcn(1,nc) = jp
                       icnpn(jp) = jp
                       ictpn(jp) = nc
                       npct(nc) = 1
                    end if
                    lactive(jp) = .true.
                 else if (iptcl == iptpn(ip)) then
                    lactive(jp) = .false.
                 end if
                 if (iptcl == iptpn(jp)) then
                    if (ictpn(jp) == 0) then
                       nc = nc + 1
                       nct = nct + 1
                       ipnsegcn(1,nc) = ip
                       ictcn(nc) = nc
                       icnpn(ip) = ip
                       ictpn(ip) = nc
                       npct(nc) = 1
                    end if
                    lactive(ip) = .true.
                 else if (iptcl == iptpn(jp)) then
                    lactive(ip) = .false.
                 end if
                 if (nbondcl(ip) == maxnbondcl(iptpn(ip))) lactive(ip) = .false.
                 if (nbondcl(jp) == maxnbondcl(iptpn(jp))) lactive(jp) = .false.
                 ! if both particles ar chain segments
              else if (lactive(ip) .and. lactive(jp)) then
                 lactive(ip) = .false.
                 lactive(jp) = .false.
                 if (bondnn(1,ip) == 0 .or. bondnn(1,jp) == 0) then
                    if (bondnn(1,ip) == 0 .and. bondnn(1,jp) == 0) then
                       icttemp = ictpn(jp)
                       nc  = nc +1
                       nct = nct + 1
                       ictpn(jp) = nc
                       icnpn(jp) = nc
                        ipnsegcn(1,icnpn(jp)) = jp
                        icttemp = ictpn(jp)
                        bondnn(1,ip) = jp
                        bondnn(1,jp) = ip
                        ipnsegcn(2,icnpn(jp)) = ip
                        icnpn(ip) = icnpn(jp)
                        ictpn(ip) = ictpn(jp)
                        npct(ictpn(jp)) = 2
                        npct(icttemp) = 0
                     else if ( bondnn(1,jp) == 0) then
                        icttemp = ictpn(jp)
                        bondnn(2,ip) = jp
                        bondnn(1,jp) = ip
                        ipnsegcn(npct(ictpn(ip)+1),icnpn(ip)) = jp
                        icnpn(jp) = icnpn(ip)
                        ictpn(jp) = ictpn(ip)
                        npct(ictpn(jp)) = npct(ictpn(jp)) + 1
                        npct(icttemp) = 0
                     else
                        icttemp = ictpn(ip)
                        bondnn(2,jp) = ip
                        ipnsegcn(npct(ictpn(jp)+1),icnpn(jp)) = ip
                        icnpn(ip) = icnpn(jp)
                        ictpn(ip) = ictpn(jp)
                        npct(ictpn(ip)) = npct(ictpn(ip)) + 1
                        npct(icttemp) = 0
                     end if

                    !bondnn(1,ip) = jp
                    !ipnsegcn(npct(ictpn(ip))+1,icnpn(ip)) = jp
                    !icnpn(jp) = icnpn(ip)
                    !ictpn(jp) = ictpn(ip)
                    !npct(ictpn(jp)) = npct(ictpn(jp)) + 1
                  else
                     bondnn(2,ip) = jp
                     bondnn(2,jp) = ip
                     ipnsegcn(npct(ictpn(ip))+1,icnpn(ip)) = jp
                     icttemp = ictpn(jp)
                     !do iseg = 1, npct(ictpn(jp))
                     !   qp = ipnsegcn(iseg, icnpn(jp))
                     !   ipnsegcn(npct(ictpn(ip) + 1 + npct(ictpn(qp)) - iseg), icnpn(ip)) = qp
                     !   icnpn(qp) = icnpn(ip)
                     !   ictpn(qp) = ictpn(ip)
                     !   npct(ictpn(ip)) = npct(ictpn(ip)) + 1
                     !end do
                     !npct(icttemp) = 0
                  end if
                     !do qp = 1, np
                       ! if (ictpn(qp) <= icttemp)  cycle
                       ! icnpn(qp) = icnpn(qp) - 1
                       ! npcttemp = npct(ictpn(qp))
                       ! ictpn(qp) = ictpn(qp) - 1
                       ! npct(ictpn(qp)) = npcttemp
                     !end do
                     !nc = nc - 1
                     !nct = nct - 1

              else
                 if (bondnn(1,ip) /= 0) then
                    bondnn(2,ip) = jp
                    ipnsegcn(npct(ictpn(ip)+1),icnpn(ip)) = jp
                    icnpn(jp) = icnpn(ip)
                    ictpn(jp) = ictpn(ip)
                    nptct(ictpn(jp)) = nptct(ictpn(jp)) + 1
                 else if (lactive(ip)) then
                    bondnn(1,ip) = jp
                    ipnsegcn(npct(ictpn(ip))+1,icnpn(ip)) = jp
                    icnpn(jp) = icnpn(ip)
                    ictpn(jp) = ictpn(ip)
                    npct(ictpn(jp)) = npct(ictpn(jp)) + 1
                 end if
                 if (bondnn(1,jp) /= 0) then
                    bondnn(2,jp) = ip
                    ipnsegcn(npct(ictpn(jp))+1,icnpn(jp)) = ip
                    icnpn(ip) = icnpn(jp)
                    ictpn(ip) = ictpn(jp)
                    npct(ictpn(ip)) = npct(ictpn(ip)) + 1
                 else if (lactive(jp)) then
                    bondnn(1,jp) = ip
                    ipnsegcn(npct(ictpn(jp))+1,icnpn(jp)) = ip
                    icnpn(ip) = icnpn(jp)
                    ictpn(ip) = ictpn(jp)
                    npct(ictpn(ip)) = npct(ictpn(ip)) + 1
                 end if
                  if (bondnn(1,ip) /= 0 .and. bondnn(2,ip) /= 0 ) then
                     lactive(ip) = .false.
                  else
                     lactive(ip) = .true.
                  end if
                  if (bondnn(1,jp) /= 0 .and. bondnn(2,jp) /= 0 ) then
                     lactive(jp) = .false.
                  else
                     lactive(jp) = .true.
                  end if
              end if
            end do
         end do
         deallocate(angle)
         deallocate(bond)
         allocate(bond(nct))
         bond = bond_var(1.0, 2, 6.0)
         allocate(angle(nct))
         angle = bond_var(0.0, 2, 0.0)


      end subroutine MakeSynthesis

      subroutine MakeSynthesis2

         use MolModule
         use Random_Module
         implicit none
         integer(4) :: ip, jp, ipt, iseg, qp, icttemp, npcttemp, i, isegqp, icntemp
         real(8) :: dx, dy, dz, rdist, Random2
         logical :: llinked

         do ip = 1, np
            if (.not. lactive(ip)) cycle
            do jp = 1, np
               if (.not. lactive(ip) ) exit
               dx = ro(1,ip) - ro(1,jp)
               dy = ro(2,ip) - ro(2,jp)
               dz = ro(3,ip) - ro(3,jp)
               call PBCr2(dx, dy, dz, rdist)
               if (rdist > 50.0) cycle   !(1.2*bondeq**2)
               if (bondnn(2,jp) /= 0) cycle
               if (ictpn(1) /= 1) print *, "error11!!! ", ictpn(1), ip, jp
               if (ictpn(jp) == ictpn(ip)) cycle
               if (ictpn(1) /= 1) print *, "error10!!! ", ictpn(1), ip, jp
               npct(ictpn(ip)) = npct(ictpn(ip)) +1
               !icttemp = ictpn(jp)
               !icntemp = icnpn(jp)
               if (ictpn(1) /= 1) print *, "error9!!! ", ictpn(1), ip, jp
               ictpn(jp) = ictpn(ip)
               if (ictpn(1) /= 1) print *, "error8!!! ", ictpn(jp), ip, jp
               ipnsegcn(npct(ictpn(ip)), icnpn(ip)) = jp
               if (ictpn(1) /= 1) print *, "error7a!!! ", ictpn(99), ip, jp
               if (ictpn(1) /= 1) print *, "error7a!!! ", icnpn(1), icnpn(99)
               if (ictpn(1) /= 1) print *, "error7a!!! ", ictpn(1), ip, jp
               if (ictpn(1) /= 1) print *, "error7a!!! ", npct(1), nc
               if (ictpn(1) /= 1) print *, "error7a!!! ", npct(ictpn(1))
               icnpn(jp) = icnpn(ip)
               if (ictpn(1) /= 1) print *, "error7!!! ", ictpn(1), ip, jp
               if (bondnn(1,ip) == 0) then
                  bondnn(1,ip) = jp
                  if (ictpn(1) /= 1) print *, "error6!!! ", ictpn(1)
               else
                  bondnn(2,ip) = jp
                  if (ictpn(1) /= 1) print *, "error5!!! ", ictpn(1)
               end if
               if (bondnn(1,jp) == 0) then
                  bondnn(1,jp) = ip
                  if (ictpn(1) /= 1) print *, "error4!!! ", ictpn(1)
               else
                  bondnn(2,jp) = ip
               if (ictpn(1) /= 1) print *, "error3!!! ", ictpn(1)
               end if
               if (lactive(jp)) then
                  if (bondnn(1,jp) /= ip .and. bondnn(1,jp) /= 0) then
                     !icttemp = ictcn(bondnn(1,jp))
                     !icntemp = icnpn(bondnn(1,jp))
                     npct(icttemp) = npct(icttemp) - 1
                     do iseg = 1, npct(icttemp)
                        if (ipnsegcn(2,icntemp) == bondnn(1,jp)) then
                           ipnsegcn(npct(ictpn(jp)) + 1, icntemp) = ipnsegcn( iseg + 1, icntemp)
                        else
                           ipnsegcn(npct(ictpn(jp)) + 1, icntemp) = ipnsegcn( npct(icttemp)+1-iseg, icntemp)
                        end if
                        ictpn(ipnsegcn(npct(ictpn(jp)) + 1, icntemp)) = ictpn(jp)
                        icnpn(ipnsegcn(npct(ictpn(jp)) + 1, icntemp)) = icnpn(jp)
                     end do
                  end if
                  npct(icttemp) = 0
                  lactive(ip) = .false.
                  lactive(jp) = .false.
                  if (ictpn(1) /= 1) print *, "error2!!! ", ictpn(1)
               else
                  lactive(ip) = .false.
                  lactive(jp) = .true.
                  if (ictpn(1) /= 1) print *, "error1!!! ", ictpn(1)
                  !do iseg = 1, npct(1)
                  !   print *, "chain segment: ", ipnsegcn(iseg,1)
                  !end do
               end if
            end do
         end do
         if (ictpn(1) /= 1) print *, "error!!! ", ictpn(1)
         !print *, "nc: ", nc
         !print *, npct(1)
      end subroutine MakeSynthesis2

      subroutine MakeSynthesis3

         use MolModule
         use Random_Module
         implicit none
         integer(4) :: ip, jp, ipt, iseg, qp, icttemp, npcttemp, i, isegqp, icntemp
         real(8) :: dx, dy, dz, rdist, Random2
         logical :: llinked

         do ip = 1, np
            if (.not. lactive(ip)) cycle
            do jp = 1, np
               if (ip == jp) cycle
               if (.not. lactive(ip) ) exit
               if (ictpn(ip) /= 0 .and. ictpn(jp) == ictpn(ip)) cycle
               if (bondnn(2,jp) /= 0) cycle
               !if (bondnn(1,jp) == ip) cycle
               if (bondnn(1,jp) /= 0 .and. nbondcl(jp) > 0) cycle
               if (bondnn(1,ip) /= 0 .and. nbondcl(ip) > 0) cycle
               if (bondcl(1,ip) == jp) cycle
               if (bondcl(1,jp) == ip) cycle
               dx = ro(1,ip) - ro(1,jp)
               dy = ro(2,ip) - ro(2,jp)
               dz = ro(3,ip) - ro(3,jp)
               call PBCr2(dx, dy, dz, rdist)
               if (rdist > 40.0) cycle   !(1.2*bondeq**2)
               if (iptcl == iptpn(ip) .and. iptcl == iptpn(jp)) cycle
               ! one particle is crosslinker
               if (iptcl == iptpn(ip) .or. iptcl == iptpn(jp)) then
                 if (nbondcl(ip) == maxnbondcl(iptpn(ip))) cycle
                 if (nbondcl(jp) == maxnbondcl(iptpn(jp))) cycle
                 nbondcl(ip) = nbondcl(ip) + 1
                 nbondcl(jp) = nbondcl(jp) + 1
                 bondcl(nbondcl(ip),ip) = jp
                 bondcl(nbondcl(jp),jp) = ip
                 if (iptcl == iptpn(ip)) then
                    if (ictpn(jp) == 0) then
                       nc = nc + 1
                       nct = nct + 1
                       ictcn(nc) = nc
                       ipnsegcn(1,nc) = jp
                       icnpn(jp) = jp
                       ictpn(jp) = nc
                       npct(nc) = 1
                       lactive(jp) = .true.
                    else
                       lactive(jp) = .false.
                    end if
                    lactive(ip) = .true.
                 end if
                 if (iptcl == iptpn(jp)) then
                    if (ictpn(ip) == 0) then
                       nc = nc + 1
                       nct = nct + 1
                       ipnsegcn(1,nc) = ip
                       ictcn(nc) = nc
                       icnpn(ip) = ip
                       ictpn(ip) = nc
                       npct(nc) = 1
                       lactive(ip) = .true.
                    else
                       lactive(ip) = .false.
                    end if
                    lactive(jp) = .true.
                 end if
                 if (nbondcl(ip) == maxnbondcl(iptpn(ip))) lactive(ip) = .false.
                 if (nbondcl(jp) == maxnbondcl(iptpn(jp))) lactive(jp) = .false.
                 ! both particles are chains
               else
                 if (bondcl(1,ip) /= 0 .and. bondcl(1,ip) == bondcl(1,jp)) cycle
                  if ( ictpn(jp) /= 0) cycle
                  npct(ictpn(ip)) = npct(ictpn(ip)) + 1
                  !icttemp = ictpn(jp)
                  !icntemp = icnpn(jp)
                  ictpn(jp) = ictpn(ip)
                  ipnsegcn(npct(ictpn(ip)), icnpn(ip)) = jp
                  icnpn(jp) = icnpn(ip)
                  if (bondnn(1,ip) == 0) then
                     bondnn(1,ip) = jp
                  else
                     bondnn(2,ip) = jp
                  end if
                  if (bondnn(1,jp) == 0) then
                     bondnn(1,jp) = ip
                  else
                     bondnn(2,jp) = ip
                  end if
                  if (lactive(jp)) then
                     lactive(ip) = .false.
                     lactive(jp) = .false.
                  else
                     lactive(ip) = .false.
                     lactive(jp) = .true.
                  end if
               end if
            end do
         end do
         !print *, "nc: ", nc
         !print *, npct(1)
      end subroutine MakeSynthesis3

      subroutine MakeSynthesis4

         use MolModule
         use Random_Module
         implicit none
         integer(4) :: ip, jp, ipt, iseg, qp, icttemp, npcttemp, i, isegqp, icntemp
         real(8) :: dx, dy, dz, rdist, Random2
         logical :: llinked

         do ip = 1, np
            if (.not. lactive(ip)) cycle
            do jp = 1, np
               if (ip == jp) cycle
               if (.not. lactive(ip) ) exit
               if (ictpn(jp) == ictpn(ip)) cycle
               !if (bondnn(2,jp) /= 0) cycle
               !if (bondnn(1,jp) /= 0) cycle
               !if (bondnn(1,jp) /= 0 .and. nbondcl(jp) > 0) cycle
               !if (bondnn(1,ip) /= 0 .and. nbondcl(ip) > 0) cycle
               !if (bondcl(1,ip) == jp) cycle
               !if (bondcl(1,jp) == ip) cycle
               dx = ro(1,ip) - ro(1,jp)
               dy = ro(2,ip) - ro(2,jp)
               dz = ro(3,ip) - ro(3,jp)
               call PBCr2(dx, dy, dz, rdist)
               if (rdist > 40.0) cycle   !(1.2*bondeq**2)
               !if (bondnn(1,jp) == ip) cycle
               !if (bondnn(2,ip) == jp) cycle
               npct(ictpn(ip)) = npct(ictpn(ip)) + 1
               ictpn(jp) = ictpn(ip)
               ipnsegcn(npct(ictpn(ip)), icnpn(ip)) = jp
               icnpn(jp) = icnpn(ip)
               if (bondnn(1,ip) == 0) then
                  bondnn(1,ip) = jp
               else
                  bondnn(2,ip) = jp
               end if
               if (bondnn(1,jp) == 0) then
                  bondnn(1,jp) = ip
               else
                  bondnn(2,jp) = ip
               end if
               if (lactive(jp)) then
                  lactive(ip) = .false.
                  lactive(jp) = .false.
               else
                  lactive(ip) = .false.
                  if (bondnn(2,jp) /= 0) then
                     lactive(jp) = .false.
                  else
                     lactive(jp) = .true.
                  end if
               end if
            end do
         end do
      end subroutine MakeSynthesis4

end module NetworkSynModule
