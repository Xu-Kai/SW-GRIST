
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: Diagnose module
! Revision history:
!----------------------------------------------------------------------------

 module swe_diagnose_module

  use grist_lib
  use grist_constants,     only: gravity, i4, r8, rearth
  use grist_domain_types
  use grist_data_types,    only: scalar_1d_field
  use grist_gcd_module,    only: calc_vorticity_at_dual_cell, &
                                 divergence_operator
  use grist_recon_module,  only: calc_scalar_at_dual_cell
  use grist_flux_operators,only: calc_normal_flux_hx_edge
  use grist_util_module,   only: write_string
  use grist_nml_module,    only: outdir, testcase, time_scheme, &
                                 advection_scheme, conserve_scheme,&
                                 pv_order, swe_timestep
  use swe_vars_module,    only:  scalar_normal_velocity_at_edge,&
                                 scalar_height_at_prime_cell   ,&
                                 scalar_topo_at_prime_cell     ,&
                                 scalar_normal_flux_at_edge
  implicit none

   private

   public   :: diagnose_total_energy,         &
               diagnose_hkinetic_energy,         &
               diagnose_potential_enstropy,   &
               diagnose_divergence_vorticity

   character(len=128)     :: c_glevel
   character(len=128)     :: c_testcase
   character(len=128)     :: c_adv
   character(len=128)     :: c_pv
   character(len=32)     :: procs

  contains

!================================================
!             Diagnose Total Energy  
!================================================

  subroutine diagnose_total_energy(mesh)

! io
    type(global_domain),  intent(in)  :: mesh
! local
    integer(i4)                        :: ie,ierr
    integer(i4)                        :: iv
    real(r8)                           :: scalar,total_energy_all1
    real(r8)                           :: total_energy,total_energy_all

!================================================
!             sum energy at edges
!================================================

     call calc_normal_flux_hx_edge(mesh,scalar_normal_velocity_at_edge%f,&  ! in
                                        scalar_normal_velocity_at_edge%f,&  ! dum
                                        scalar_height_at_prime_cell%f   ,&  ! in
                                        scalar_normal_flux_at_edge%f    ,&  ! out
                                        advection_scheme                ,&
                                        swe_timestep )

     scalar_normal_flux_at_edge%f = scalar_normal_flux_at_edge%f*scalar_normal_velocity_at_edge%f
!
! sum energy at edges 
!
     total_energy = 0._r8
     !do ie = 1, mesh%ne
     !   total_energy = total_energy+mesh%edp(ie)%leng*mesh%edt(ie)%leng*(rearth**2)*&
     !                               scalar_normal_flux_at_edge%f(ie)*&
     !                               scalar_normal_velocity_at_edge%f(ie)*0.5_r8
     !end do
     do ie = 1, mesh%ne
        if(mesh%edt(ie)%v(1)>mesh%nv_compute.or.mesh%edt(ie)%v(2)>mesh%nv_compute)then
                total_energy = total_energy+mesh%edp(ie)%leng*mesh%edt(ie)%leng*(rearth**2)*&
                                        scalar_normal_flux_at_edge%f(ie)*&
                                        scalar_normal_velocity_at_edge%f(ie)*0.25_r8
        else
                total_energy = total_energy+mesh%edp(ie)%leng*mesh%edt(ie)%leng*(rearth**2)*&
                                        scalar_normal_flux_at_edge%f(ie)*&
                                        scalar_normal_velocity_at_edge%f(ie)*0.5_r8
        endif
     end do
!
! sum energy at cells
!
     do iv = 1, mesh%nv
        scalar       = scalar_height_at_prime_cell%f(iv)
        total_energy = total_energy+gravity*scalar*&
                      (scalar/2._r8+scalar_topo_at_prime_cell%f(iv))*((rearth**2)*mesh%plg(iv)%areag)
     end do
     call reduce(total_energy, total_energy_all, 'sum')

     call write_string(mesh%glevel,c_glevel)
     call write_string(testcase   ,c_testcase)
     call write_string(advection_scheme,c_adv)
     call write_string(pv_order(3),c_pv)

     if(mpi_rank()==0)then
        open(1,file=trim(outdir)//&
                "GRIST-TC"//trim(c_testcase)//&
                "-GLEVEL"//trim(c_glevel)//&
                "-"//trim(time_scheme)//&
                "-PV"//trim(c_pv)//&
                "-ADV"//trim(c_adv)//&
                "-"//trim(conserve_scheme)//&
                "-TOTAL-ENERGY.txt",access='append')
        write(1,*) total_energy_all
        close(1)
     end if
     
    return
  end subroutine diagnose_total_energy

!================================================
! Calculate H*kinetic energy at prime cell
!================================================

  subroutine diagnose_hkinetic_energy(mesh)
! io
   type(global_domain), intent(in)    :: mesh
! local
   integer                             :: iv,ie,ierr
   real(r8)                            :: length_of_voronoi_edge
   real(r8)                            :: length_of_triangle_edge
   real(r8)                            :: tmp_sum
   real(r8)                            :: wind
   real(r8)                            :: hkinetic_energy,hkinetic_energy_all

   hkinetic_energy = 0._r8

   do iv=1, mesh%nv
      tmp_sum  = 0._r8
      do ie = 1, mesh%vtx(iv)%nnb
           length_of_voronoi_edge  = rearth*mesh%edp(mesh%vtx(iv)%ed(ie))%leng
           length_of_triangle_edge = rearth*mesh%edt(mesh%vtx(iv)%ed(ie))%leng
           wind                    = scalar_normal_velocity_at_edge%f(mesh%vtx(iv)%ed(ie))
           tmp_sum                 = tmp_sum+0.25_r8*length_of_voronoi_edge*length_of_triangle_edge*(wind**2)
      end do
           hkinetic_energy    = hkinetic_energy+scalar_height_at_prime_cell%f(iv)*tmp_sum
   end do
   call reduce(hkinetic_energy, hkinetic_energy_all, 'sum')

     call write_string(mesh%glevel,c_glevel)
     call write_string(testcase   ,c_testcase)
     call write_string(advection_scheme,c_adv)
     call write_string(pv_order(3),c_pv)

     if(mpi_rank()==0)then
       open(1,file=trim(outdir)//&
                 "GRIST-TC"//trim(c_testcase)//&
                 "-GLEVEL"//trim(c_glevel)//&
                 "-"//trim(time_scheme)//&
                 "-PV"//trim(c_pv)//&
                 "-ADV"//trim(c_adv)//&
                 "-"//trim(conserve_scheme)//&
                 "-HKENERGY.txt",access='append')
       write(1,*) hkinetic_energy_all
       close(1)
     end if

   return
  end subroutine diagnose_hkinetic_energy

!================================================
!           Diagnose Potential Enstropy 
!================================================

  subroutine diagnose_potential_enstropy(mesh)
! io
    type(global_domain),  intent(inout) :: mesh
! local
    integer(i4)                          :: it,ierr
    logical                              :: f1,f2,f3
    real(r8)                             :: scalar,potential_enstropy_all
    real(r8)                             :: potential_enstropy
    type(scalar_1d_field)                :: scalar_height_at_dual_cell
    type(scalar_1d_field)                :: scalar_absolute_vorticity_at_dual_cell
    type(scalar_1d_field)                :: scalar_relative_vorticity_at_dual_cell

     allocate(scalar_height_at_dual_cell%f(1:mesh%nt_full))
     allocate(scalar_absolute_vorticity_at_dual_cell%f(1:mesh%nt_full))
     allocate(scalar_relative_vorticity_at_dual_cell%f(1:mesh%nt_full))

     call calc_scalar_at_dual_cell(mesh,scalar_height_at_prime_cell%f,&
                                        scalar_height_at_dual_cell%f)

     call calc_vorticity_at_dual_cell(mesh,scalar_normal_velocity_at_edge%f         ,&
                                           scalar_absolute_vorticity_at_dual_cell%f ,&
                                           scalar_relative_vorticity_at_dual_cell%f )
     potential_enstropy = 0._r8

     !do it = 1, mesh%nt_iner
     !    potential_enstropy = potential_enstropy+&
     !                        (rearth**2)*mesh%tri(it)%areag*&
     !                        (scalar_absolute_vorticity_at_dual_cell%f(it)**2)/&
     !                         scalar_height_at_dual_cell%f(it)
     !end do
      
     do it = 1, mesh%nt
        f1=mesh%tri(it)%v(1)<=mesh%nv_compute
        f2=mesh%tri(it)%v(2)<=mesh%nv_compute
        f3=mesh%tri(it)%v(3)<=mesh%nv_compute
        if(((f1).and.(.not.f2).and.(.not.f3)).or.((.not.f1).and.&
        (f2).and.(.not.f3)).or.((.not.f1).and.(.not.f2).and.(f3)))then
         potential_enstropy = potential_enstropy+&
                             (rearth**2)*mesh%tri(it)%areag*&
                             (scalar_absolute_vorticity_at_dual_cell%f(it)**2)/&
                              scalar_height_at_dual_cell%f(it)*(1/3._r8)
        elseif(((f1).and.(f2).and.(.not.f3)).or.((f1).and.(.not.f2).and.&
        (f3)).or.((.not.f1).and.(f2).and.(f3))) then
         potential_enstropy = potential_enstropy+&
                             (rearth**2)*mesh%tri(it)%areag*&
                             (scalar_absolute_vorticity_at_dual_cell%f(it)**2)/&
                              scalar_height_at_dual_cell%f(it)*(2/3._r8)
        else
         potential_enstropy = potential_enstropy+&
                             (rearth**2)*mesh%tri(it)%areag*&
                             (scalar_absolute_vorticity_at_dual_cell%f(it)**2)/&
                              scalar_height_at_dual_cell%f(it)
        endif
     end do
     call reduce(potential_enstropy, potential_enstropy_all, 'sum') 

     deallocate(scalar_height_at_dual_cell%f)
     deallocate(scalar_absolute_vorticity_at_dual_cell%f)
     deallocate(scalar_relative_vorticity_at_dual_cell%f)

     call write_string(mesh%glevel,c_glevel)
     call write_string(testcase   ,c_testcase)
     call write_string(advection_scheme,c_adv)
     call write_string(pv_order(3),c_pv)

     if(mpi_rank()==0)then
       open(1,file=trim(outdir)//&
                 "GRIST-TC"//trim(c_testcase)//&
                 "-GLEVEL"//trim(c_glevel)//&
                 "-"//trim(time_scheme)//&
                 "-PV"//trim(c_pv)//&
                 "-ADV"//trim(c_adv)//&
                 "-"//trim(conserve_scheme)//&
                 "-POTENTIAL-ENSTROPHY.txt",access='append')
       write(1,*) potential_enstropy_all
       close(1)
     end if

     return
   end subroutine diagnose_potential_enstropy

!================================================
!    Diagnose vorticity and divergence
!================================================

   subroutine diagnose_divergence_vorticity(mesh,scalar_normal_velocity_at_edge,        &
                                                 scalar_height_at_prime_cell,           &
                                                 scalar_divergence_at_prime_cell,       &
                                                 scalar_absolute_vorticity_at_dual_cell,&
                                                 scalar_relative_vorticity_at_dual_cell,&
                                                 scalar_potential_vorticity_at_dual_cell)
! io
     type(global_domain),      intent(in)    :: mesh
     type(scalar_1d_field),    intent(in)    :: scalar_normal_velocity_at_edge
     type(scalar_1d_field),    intent(in)    :: scalar_height_at_prime_cell
     type(scalar_1d_field),    intent(inout) :: scalar_divergence_at_prime_cell
     type(scalar_1d_field),    intent(inout) :: scalar_absolute_vorticity_at_dual_cell
     type(scalar_1d_field),    intent(inout) :: scalar_relative_vorticity_at_dual_cell
     type(scalar_1d_field),    intent(inout) :: scalar_potential_vorticity_at_dual_cell
! local
     type(scalar_1d_field)                   :: scalar_height_at_dual_cell
     type(scalar_1d_field)                   :: scalar_height_unity         ! set to one
     integer(i4)                             :: it

     allocate(scalar_height_at_dual_cell%f(1:mesh%nt_full))

     call divergence_operator(mesh, scalar_normal_velocity_at_edge%f, &
                                    scalar_divergence_at_prime_cell%f )
   
     call calc_scalar_at_dual_cell(mesh,scalar_height_at_prime_cell%f,&
                                        scalar_height_at_dual_cell%f)

     call calc_vorticity_at_dual_cell(mesh,scalar_normal_velocity_at_edge%f        ,&
                                           scalar_absolute_vorticity_at_dual_cell%f,&
                                           scalar_relative_vorticity_at_dual_cell%f )

     do it = 1, mesh%nt
        scalar_potential_vorticity_at_dual_cell%f(it) = scalar_absolute_vorticity_at_dual_cell%f(it)/&
                                                        scalar_height_at_dual_cell%f(it)
     end do

     deallocate(scalar_height_at_dual_cell%f)

    return
   end subroutine diagnose_divergence_vorticity

  end module swe_diagnose_module
