!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description:
!     Definitions of discrete elements. Some elements are saved by 
! grid_generator (gg) or mpi_scvt, marked by "saved by gg". Other are evaluated 
! in the initilization process. Note that the dual cell is currently limited 
! to a triangle, that is, a Voronoi-Delaunay association. It can be generalized 
! to a polygonal cell. But given the current Staggered-FV numerics, there is no
! reason to deviate from an SCVT-optimized mesh, which is perhaps the best choice
! Note: all these element are deleted in 3D model runtime.
!
! Revision history:
!       1. For cubed sphere grid, some mesh weight may be obsolete as they 
! are specific to icosahedral grid; some weight may require rewrite
!----------------------------------------------------------------------------

 module grist_element_types_cubed

  !use grist_constants, only: i4, r8
  use grist_constants_dbl, only: i4, r8

  implicit none

    public :: vector               , &
              node_structure       , & ! original generated node (point)
              edge_structure       , & ! edge (line)
              dual_structure       , & ! dual (cell/tr), hard coded to triangle now
              prim_structure       , & ! prime(cell/hx), arbitrarily-structured polygon
              node_structure_global, &
              edge_structure_global, &
              dual_structure_global

    integer(i4), parameter :: maxnnb  = 4  ! normally, only 7-sided polygons are found
    integer(i4), parameter :: maxnnb2 = 8  ! maxnnb*2
    integer(i4), parameter :: nTwo    = 2
    integer(i4), parameter :: nThree  = 3
    integer(i4), parameter :: nDual   = 4  ! only account for edge neighbors

  type vector
     real(r8), dimension(1:nThree) :: v
  end type vector

  type point_type
     real(r8), dimension(1:nThree) :: p           ! cartesian coordinates in vector form
     real(r8)                      :: lat, lon    ! spherical coordinates of node in lat [-pi/2, pi/2], lon [-pi, pi]
     real(r8), dimension(1:nThree) :: north       ! a unit vector pointing to the local nort direction
     real(r8), dimension(1:nThree) :: east        ! a unit vector pointing to the local east direction
  end type point_type

!------------------------------------------------------------------
! Triangle vertice/Voronoi cell node structure
! Node is the vertice of triangle, i.e., original SCVT/icosahedral
! generators
!------------------------------------------------------------------

  type node_structure

! below 7 saved by gg (5 vtx_xx)
     real(r8)                 :: p(nThree)   ! cartesian coordinates in vector form
     real(r8)                 :: lon, lat    ! spherical coordinates of node in lat [-pi/2, pi/2], lon [-pi, pi]
     integer(i4)              :: nnb         ! connectivity: number of neighbour nodes (conected through edges)
     integer(i4), allocatable :: nb(:)       ! connectivity: index of neighbour vertices; counterclockwise sequences of neighboring nodes
     integer(i4), allocatable :: ed(:)       ! index of edges that has this point as vertice; size: 1:nnb
     integer(i4), allocatable :: tr(:)       ! index of triangles/corner that has this point as vertice; size: 1:nnb
!
! local orthogonal basis (lob)'s nx and ny, used by FFSL-style flux
!
     real (r8), allocatable   :: lob_nx(:,:) ! each lob-X connects this cell with one of its neighbor
     real (r8), allocatable   :: lob_ny(:,:) ! each lob-Y is 90 counterclockwise of lob-X

  end type node_structure

!------------------------------------------------------------------
! an edge may belong to a dual or a prim cell
!------------------------------------------------------------------

  type edge_structure

! below 6 saved by gg (5 edt/edp_, +ltln is provided for plotting, not use here, similar for below)
     type(point_type)             :: c       ! crossing point between triangle and voronoi edge, saved by gg
     integer (i4), dimension(1:nTwo) :: v       ! global index of end points of this edge, saved by gg
     real(r8)                     :: lenp    ! projected length- Euclidian distance in R3 considers the unit sphere, saved by gg
     real(r8)                     :: leng    ! geodesical length considering the unit sphere, saved by gg
!
! the tangen unit vector of edge, placed at the edge point, tangent to the sphere,
! for edt: edt's nr x k
! for edp: edp's nr x k
! saved by gg
!
     real (r8), dimension(1:nThree)  :: tg

!
! the normal unit vector of edge, placed at the edge point , tangent to the sphere, 
! for edt: cross product of two ends (v1 X v2)
! for edp: = edt's tg， i.e., v(1)->v(2)
! saved by gg
!
     real (r8), dimension(1:nThree)  :: nr
!
! only for edp, the triangular points that constitute a pair of tr shared by edt, non edt's end point
!
     integer (i4), dimension(1:2)    :: vtx
!
! indicator function to check plg's nr direction relative to tr's tg direction 
! dot_product(mesh%edp(ie)%nr,mesh%edt(ie)%tg)
! If they are the same direction,   edpNr_edtTg= 1
! If they are in opposite direction,edpNr_edtTg=-1
! given the present config, this is always 1, but we explicity code it in case 
! that the precondition changes
!
     integer(i4)  :: edpNr_edtTg
!
! indicator function to check plg's tg direction relative; to tr's nr direction, only defined 
! by grist; =dot_product(mesh%edp(ie)%tg,mesh%edt(ie)%nr)
! If they are the same direction,   edpTg_edtNr= 1
! If they are in opposite direction,edpTg_edtNr=-1
!
     integer(i4)  :: edpTg_edtNr
!
! indicator to check whether edp's v1->v2 direction coincide with the tg direction of edp
! dot_product(edp(ie)%v2-v1,mesh%edp(ie)%tg)
! If they are the same direction,    edp12_edp_tg = 1
! If they are in opposite direction, edp12_edp_tg =-1
!
     integer(i4)  :: edp12_edp_tg
!
! indicator to check whether edp's v1->v2 direction coincide with the nr direction of edt
! dot_product(edp(ie)%v2-v1,mesh%edp(ie)%nr)
! If they are the same direction,    edp12_edt_nr = 1
! If they are in opposite direction, edp12_edt_nr =-1
!
     integer(i4)  :: edp12_edt_nr

!
! local index of edge in the plg or tri used by the coriolis term calculation
!
     integer(i4)  :: local_index_prime_cell(nTwo)
     integer(i4)  :: local_index_dual_cell(nTwo)
!
! Local weights used by the polynomial reconstruction; they are stored by infering 
! those stored in the face-element (plg/tri), but indexed by edge, each end point of an edge 
! has such a weight, the same local index as end 'v' 
! avoid allocatable matrix as suggested
!

! for 2nd polynomial
     real(r8)                     :: v_weight_1st_prime_cell(nTwo,maxnnb)   ! assume maxnnb=8
     real(r8)                     :: v_weight_prime_cell(nTwo,maxnnb)
! for 4th polynomial
     real(r8)                     :: v_weight_4th_prime_cell(nTwo,23)
     real(r8)                     :: v_weight_4th_2nd_prime_cell(nTwo,23)
! mid point location of edge
     real(r8)                     :: v_mid_point_lob_prime_cell(nTwo,nTwo) ! xy in v1v0's lob
! for 2nd polynomial
     real(r8)                     :: v_weight_dual_cell(nTwo,maxnnb2)
! for 4th polynomial, experimental (not used since 2018)
     real(r8)                     :: v_weight_4th_dual_cell(nTwo,26)
! mid point location of edge 
     real(r8)                     :: v_mid_point_lob_dual_cell(nTwo,nTwo)  ! xy in v1v0's lob
! for runtime nct code
     integer (i4)                 :: edge_on_edge(maxnnb2)  ! 16 as max, edge index
     real (r8)                    :: trsk_on_edge(maxnnb2)  ! trsk weight
     real (r8)                    :: edpl_on_edge(maxnnb2)  ! edp's length
     integer(i4)                  :: nedge             ! how many edge on edge
! for runtime O58
     integer(i4)                  :: my_edge_on_edge(nTwo,maxnnb) ! dim1/2: global edge index on v1/2 that share this edge, begin from this edge, ccw
     integer(i4)                  :: ur_cell_on_edge(nTwo,maxnnb) ! dim1/2: local cell index on one nb ed of v1/2 that share this edge, begin from this edge's nb cell, ccw
     integer(i4)                  :: my_edge_on_edge_num(nTwo)    ! actual number of second dim of above array
     integer(i4)                  :: ur_cell_on_edge_num(nTwo)
!
! for Gassmann PV&NCT method (testing implementation around 2018, deleted now), has not been used (for reference only!)
!
     integer(i4)                  :: edge_on_rhombi(4) ! global index of edge
     integer(i4)                  :: eccw_on_rhombi(4) ! check whether normal direction of edp conforms to the CCW direction relative to the target edge point
     real(r8)                     :: edtl_on_rhombi(4) ! tr edge length
     real(r8)                     :: area_on_rhombi    ! rhombi area
     integer(i4)                  :: gass_edge_on_edge(maxnnb2) ! gassmann's edge index on edge for NCT
     integer(i4)                  :: gass_flag_on_edge(maxnnb2) ! gassmann's flag index on edge for NCT

  end type edge_structure

!=================================================================
! triangle/dual cell. Note that we assume it is three-sided, 
! this assumption can be relaxed in future, but would require a
! significant change of parallel infrastructure as well
!=================================================================

  type dual_structure

! below 8 saved by gg (7 tri_xx+ ltln)
     type(point_type)                :: c     ! circumcenter point, saved by gg
     real(r8)                        :: areag ! area of geodesical triangle of unit sphere, saved by gg
     integer(i4), dimension(1:nDual) :: v     ! index of generating vertices; counter clockwise ordering, saved by gg
     integer(i4), dimension(1:nDual) :: ed    ! index of triangle edges, the ith edge is conecting v(i) and v(i+1), saved by gg
     integer(i4), dimension(1:nDual) :: nb    ! index of edge-neighbour triangles, saved by gg
     integer(i4)                     :: nnb   ! added for cubed-sphere mesh; for icos, this is always 3
!
! a binary indicator such that for each edge,
! = 1 if the direction of tangent vector CCW to cell
! =-1 else; saved by gg
!
     integer(i4), dimension(1:nDual) :: tg
!
! a binary indicator such that for each edge,
! 1 if the direction of normal vector points out of cell; -1 else
! saved by gg
!
     integer(i4), dimension(1:nDual) :: nr
!
! dual-primal cell intersection kite area, belonged by triangle
!
     real(r8), dimension(1:nDual):: kite_area
!
! B array, size is m,n,m, where m is the number of neighbouring cells, which is 3(nDual), n is 5/14,
! each plg has one B for each of its edges, and each edge has two for each of its vertices.
! also defined in ed structure.
!
     real(r8),   allocatable     :: blocal(:,:,:)
     real(r8),   allocatable     :: blocal_4th(:,:,:)
     real(r8)                    :: mid_point_lob(nThree,nThree,nTwo) ! 3 lob basis * 3 edges * 2 directions (X/Y) 
!
! for local stencil reconstruction
!
     integer(i4), allocatable    :: stencil_index_2nd(:) ! index  of 2nd-order stencil
     integer(i4), allocatable    :: stencil_index_4th(:) ! index  of 4th-order stencil
     integer(i4)                 :: stencil_number_2nd   ! number of 2nd-order stencil
     integer(i4)                 :: stencil_number_4th   ! number of 4th-order stencil
!
! local orthogonal basis (lob)'s nx and ny, only used by FFSL-style flux
!
     real (r8), allocatable   :: lob_nx(:,:) ! each lob-X connects this cell with one of its neighbor
     real (r8), allocatable   :: lob_ny(:,:) ! each lob-Y is 90 counterclockwise of lob-X

  end type dual_structure

!=================================================================
! structure of polygon/voronoi/hexagon/pentagon face, each face is
! centered at on grid vertice/node so that the vertice index shares
! with this one
!=================================================================

  type prim_structure

! below 5 saved by gg (4 plg_xx+ltln)
    !type(point_type)         :: b          ! barycenter Point, not vtx'p, saved by gg
     real(r8)                 :: areag       ! area of spherical polygon, saved by gg
     integer(i4), allocatable :: tg(:)       ! a binary indicator, 1 if the direction of tangen vector CCW to cell; -1 else; no-longer saved by gg
     integer(i4), allocatable :: nr(:)       ! a binary indicator, 1 if the direction of normal vector points out of cell; -1 else; no-longer saved by gg

!
! trsk weight for vector reconstructure (cf., Thuburn et al. 2009)
!
     real(r8), allocatable    :: trsk_weight(:,:)

!
! Dual-Primal cell intersection area, belonged by polygon
!
     real (r8), allocatable   :: kite_area(:)
!
! B array , size is m,n,stencil_number, where m is the number of neighbouring cells, 
! n is 5/14, each polygon has m=1 means that this blocal is formed by the lob defined by 
! connecting THIS cell and its inb=1 neighbor one B for each of its edges, and each 
! edge has two for each of its vertices. Also defined in the ed structure
!
     real(r8),    allocatable :: blocal(:,:,:)        ! m*n*stencil_number
     real(r8),    allocatable :: blocal_4th(:,:,:)
     real(r8),    allocatable :: mid_point_lob(:,:,:) ! m*m*2
!
! stencil
!
     integer(i4), allocatable :: stencil_index_2nd(:) ! index of 2nd-order stencil
     integer(i4), allocatable :: stencil_index_4th(:) ! index of 4th-order stencil
     integer(i4)              :: stencil_number_2nd   ! number of 2nd-order stencil
     integer(i4)              :: stencil_number_4th   ! number of 4th-order stencil

!
! Midpoint locations of sub-triangle's edges of this polygon, dim(nnb,3,3)<=>dim(i,j,k)
! 1st dim denotes number of neighbors; 2nd dim denotes 3 midpoints of triangle's 3 edges
! plus area of this sub triangle; 3rd dim denotes 3D cartesian coordinate and two 
! locations in lob_x/lob_y; For example:
! i=1,   means tr connects this cell, tr1 and tr2
! i=nnb, means tr connects this cell, trnnb and tr1
! j=1, means edge connects this cell and tr1
! j=2, means edge connects this cell and tr2
! j=3, means edge connects tr1 and tr2
! see Fig. 2a in SM10 for details
!
      real(r8), allocatable   :: sub_triangle_midp_3d(:,:,:)   ! nnb*3edges*3cartedian
      real(r8), allocatable   :: sub_triangle_midp_lob(:,:,:)  ! nnb*3edges*2(lob_x,lob_y)
      real(r8), allocatable   :: sub_triangle_area(:)

  end type prim_structure

!=================================================================
! triangle vertice/raw node differ from up, we only need a part of
! node_structure
!=================================================================

  type node_structure_global

     real(r8)                 :: lon, lat ! spherical coordinates of node in radians lat in [-pi/2, pi/2] , lon in [-pi, pi]   
     integer(i4)              :: nnb      ! number of neighbour nodes (conected through edges)
     integer(i4), allocatable :: nb(:)    ! index of neighbour vertices; counterclockwise sequences of neighboring nodes
     integer(i4), allocatable :: ed(:)    ! index of edges that has this point as vertice; size: 1:nnb
     integer(i4), allocatable :: tr(:)    ! index of triangles that has this point as vertice; size: 1:nnb

!================================================
! local orthogonal basis (lob)'s nx and ny
! only used by FFSL-style flux
!================================================

     !real (r8), allocatable   :: lob_nx(:,:) ! each lob-X connects this cell with one of its neighbor
     !real (r8), allocatable   :: lob_ny(:,:) ! each lob-Y is 90 counterclockwise of lob-X

  end type node_structure_global

!search io needs this struct
  type dual_structure_global
     integer(i4), dimension(1:nDual) :: v  ! index of generating vertices; counter clockwise ordering, saved by gg
  end type dual_structure_global

!search io needs this struct
  type edge_structure_global
     integer (i4), dimension(1:2) :: v     ! global index of end points of this edge, saved by gg 
  end type edge_structure_global

 end module grist_element_types_cubed
