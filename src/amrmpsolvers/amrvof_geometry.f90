!> AMR VOF Geometry module
!> Provides cutting tables and volume computation for native VOF geometry
module amrvof_geometry
   use precision, only: WP
   implicit none
   private

   ! Expose tables and routines
   public :: tet_map, cut_side, cut_v1, cut_v2, cut_vtet
   public :: cut_ntets, cut_nvert, cut_nntet
   public :: get_plane_dist
   public :: tet_vol, tet_sign, cut_tet_vol
   public :: flux_polyhedron_vol, cut_hex_vol
   public :: correct_flux_poly
   public :: cut_hex_polygon, hex_poly_nvert, get_hex_poly_nvert
   public :: remap_box_staggered,remap_box_collocated
   public :: flux_poly_moments
   public :: project_rk2_collocated,project_rk2_staggered
   public :: tet2flux_plic,tet2flux_recursive,tet2flux_iterative

   ! Cutting tables from mpcomp_class_noirl
   ! tet_map: maps a hex cell (8 vertices + center) to 8 tetrahedra
   integer, dimension(4,8), parameter :: tet_map = reshape([ &
      7, 4, 3, 6, 6, 3, 2, 4, 6, 2, 1, 4, 7, 8, 4, 6, &
      6, 5, 8, 4, 6, 5, 4, 1, 5, 6, 8, 9, 6, 7, 8, 9], shape(tet_map))

   ! cut_side: which side of plane each tet vertex is on (1=below, 2=above)
   integer, dimension(6,16), parameter :: cut_side = reshape([ &
      1,-1,-1,-1,-1,-1, 2, 1, 1, 1,-1,-1, 2, 1, 1, 1,-1,-1, 2, 2, 2, 1, 1, 1, &
      2, 1, 1, 1,-1,-1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1,-1,-1, &
      2, 1, 1, 1,-1,-1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1,-1,-1, &
      2, 2, 2, 1, 1, 1, 2, 2, 2, 1,-1,-1, 2, 2, 2, 1,-1,-1, 2,-1,-1,-1,-1,-1], shape(cut_side))

   ! cut_v1, cut_v2: edge endpoints for cut vertices
   integer, dimension(4,16), parameter :: cut_v1 = reshape([ &
      -1,-1,-1,-1, 1, 1, 1,-1, 2, 2, 2,-1, 1, 2, 1, 2, 3, 3, 3,-1, 1, 3, 1, 3, &
      2, 3, 2, 3, 4, 4, 4,-1, 4, 4, 4,-1, 1, 4, 1, 4, 2, 4, 2, 4, 3, 3, 3,-1, &
      3, 4, 3, 4, 2, 2, 2,-1, 1, 1, 1,-1,-1,-1,-1,-1], shape(cut_v1))

   integer, dimension(4,16), parameter :: cut_v2 = reshape([ &
      -1,-1,-1,-1, 2, 3, 4,-1, 3, 4, 1,-1, 4, 4, 3, 3, 4, 1, 2,-1, 4, 4, 2, 2, &
      4, 4, 1, 1, 1, 2, 3,-1, 1, 2, 3,-1, 3, 3, 2, 2, 3, 3, 1, 1, 4, 1, 2,-1, &
      2, 2, 1, 1, 3, 4, 1,-1, 2, 3, 4,-1,-1,-1,-1,-1], shape(cut_v2))

   ! cut_vtet: sub-tetrahedra for each cutting configuration
   integer, dimension(4,6,16), parameter :: cut_vtet = reshape([ &
      1, 2, 3, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &
      5, 7, 6, 1, 6, 2, 3, 4, 4, 2, 5, 6, 5, 6, 7, 4,-1,-1,-1,-1,-1,-1,-1,-1, &
      7, 5, 6, 2, 1, 3, 4, 6, 1, 5, 3, 6, 5, 7, 6, 1,-1,-1,-1,-1,-1,-1,-1,-1, &
      5, 8, 6, 2, 5, 7, 8, 1, 5, 1, 8, 2, 5, 6, 8, 4, 5, 8, 7, 3, 5, 8, 3, 4, &
      6, 5, 7, 3, 2, 1, 4, 6, 6, 5, 4, 2, 6, 7, 5, 2,-1,-1,-1,-1,-1,-1,-1,-1, &
      5, 6, 8, 3, 5, 8, 7, 1, 5, 8, 1, 3, 5, 8, 6, 4, 5, 7, 8, 2, 5, 8, 4, 2, &
      8, 6, 5, 3, 5, 7, 8, 2, 8, 5, 2, 3, 8, 5, 6, 4, 5, 8, 7, 1, 5, 8, 1, 4, &
      1, 2, 3, 7, 1, 2, 7, 6, 5, 7, 6, 1, 5, 6, 7, 4,-1,-1,-1,-1,-1,-1,-1,-1, &
      5, 6, 7, 4, 1, 2, 3, 6, 5, 1, 3, 6, 5, 7, 6, 3,-1,-1,-1,-1,-1,-1,-1,-1, &
      5, 8, 6, 4, 5, 7, 8, 1, 5, 8, 4, 1, 5, 6, 8, 3, 5, 8, 7, 2, 5, 8, 2, 3, &
      8, 5, 6, 4, 5, 8, 7, 2, 8, 2, 5, 4, 8, 6, 5, 3, 5, 7, 8, 1, 5, 8, 3, 1, &
      1, 4, 2, 7, 4, 1, 6, 7, 6, 7, 5, 4, 6, 5, 7, 3,-1,-1,-1,-1,-1,-1,-1,-1, &
      8, 6, 5, 4, 5, 7, 8, 3, 8, 4, 5, 3, 8, 5, 6, 2, 5, 8, 7, 1, 5, 8, 1, 2, &
      3, 4, 1, 7, 7, 6, 3, 4, 7, 6, 5, 3, 7, 5, 6, 2,-1,-1,-1,-1,-1,-1,-1,-1, &
      7, 4, 2, 3, 2, 3, 6, 7, 5, 6, 7, 2, 5, 7, 6, 1,-1,-1,-1,-1,-1,-1,-1,-1, &
      1, 2, 3, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1], shape(cut_vtet))

   ! Number of resulting tetrahedra for each config
   integer, dimension(16), parameter :: cut_ntets = [1,4,4,6,4,6,6,4,4,6,6,4,6,4,4,1]
   ! Number of cut vertices for each config
   integer, dimension(16), parameter :: cut_nvert = [0,3,3,4,3,4,4,3,3,4,4,3,4,3,3,0]
   ! Number of resulting tetrahedra on one side
   integer, dimension(16), parameter :: cut_nntet = [1,2,2,4,2,4,4,4,2,4,4,4,4,4,4,2]

   !> ============================================================================
   !> Hex polygon cutting tables (from IRL lookup_tables.h)
   !> For cutting a plane by a hexahedron to get the intersection polygon
   !> Vertex ordering (x,y,z):
   !>   0:(+,-,-), 1:(+,+,-), 2:(+,+,+), 3:(+,-,+)
   !>   4:(-,-,-), 5:(-,+,-), 6:(-,+,+), 7:(-,-,+)
   !> ============================================================================
   
   ! Number of polygon vertices for each of 256 cases
   ! -1 indicates invalid/ambiguous case
   integer, dimension(0:255), parameter :: hex_poly_nvert = [ &
      0,  3,  3,  4,  3, -1,  4,  5,  3,  4, -1,  5,  4,  5,  5,  4, &
      3,  4, -1,  5, -1, -1, -1, -1, -1,  5, -1,  6, -1, -1, -1,  5, &
      3, -1,  4,  5, -1, -1,  5,  6, -1, -1, -1, -1, -1, -1, -1,  5, &
      4,  5,  5,  4, -1, -1, -1,  5, -1, -1, -1,  5, -1, -1, -1,  4, &
      3, -1, -1, -1,  4, -1,  5, -1, -1, -1, -1, -1,  5, -1,  6,  5, &
     -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
      4, -1,  5, -1,  5, -1,  4,  5, -1, -1, -1, -1, -1, -1,  5,  4, &
      5, -1,  6,  5, -1, -1,  5,  4, -1, -1, -1, -1, -1, -1, -1,  3, &
      3, -1, -1, -1, -1, -1, -1, -1,  4,  5, -1, -1,  5,  6, -1,  5, &
      4,  5, -1, -1, -1, -1, -1, -1,  5,  4, -1,  5, -1,  5, -1,  4, &
     -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
      5,  6, -1,  5, -1, -1, -1, -1, -1,  5, -1,  4, -1, -1, -1,  3, &
      4, -1, -1, -1,  5, -1, -1, -1,  5, -1, -1, -1,  4,  5,  5,  4, &
      5, -1, -1, -1, -1, -1, -1, -1,  6,  5, -1, -1,  5,  4, -1,  3, &
      5, -1, -1, -1,  6, -1,  5, -1, -1, -1, -1, -1,  5, -1,  4,  3, &
      4,  5,  5,  4,  5, -1,  4,  3,  5,  4, -1,  3,  4,  3,  3,  0]

   ! Edge pairs for polygon vertices: hex_poly_edges(1:2, vertex, case)
   ! Each pair (v1,v2) defines an edge; interpolate to find intersection point
   ! Stored as hex_poly_v1 and hex_poly_v2 for the two endpoints
   integer, dimension(6,0:255), parameter :: hex_poly_v1 = reshape([ &
      -1,-1,-1,-1,-1,-1,  0, 3, 0,-1,-1,-1,  0, 1, 1,-1,-1,-1,  1, 3, 0, 1,-1,-1, &
       1, 2, 2,-1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 1, 2, 2,-1,-1,  2, 3, 0, 1, 2,-1, &
       2, 3, 3,-1,-1,-1,  0, 2, 3, 0,-1,-1, -1,-1,-1,-1,-1,-1,  1, 2, 3, 0, 1,-1, &
       1, 2, 3, 3,-1,-1,  0, 1, 2, 3, 0,-1,  0, 1, 2, 3, 3,-1,  0, 1, 2, 3,-1,-1, &
       4, 0, 7,-1,-1,-1,  0, 3, 7, 4,-1,-1, -1,-1,-1,-1,-1,-1,  1, 3, 7, 4, 1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1,  0, 2, 3, 7, 4,-1, -1,-1,-1,-1,-1,-1,  1, 2, 3, 7, 4, 1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  4, 1, 2, 3, 7,-1, &
       4, 5, 1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 4, 5, 1,-1,-1,  1, 3, 0, 4, 5,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 4, 5, 2, 2,-1,  2, 3, 0, 4, 5, 2, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  4, 5, 2, 3, 0,-1, &
       5, 1, 0, 7,-1,-1,  0, 3, 7, 5, 1,-1,  0, 0, 7, 5, 1,-1,  1, 3, 7, 5,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 5, 2,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 2, 3, 7, 5,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  5, 2, 3, 7,-1,-1, &
       5, 6, 2,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 5, 6, 2,-1,-1, -1,-1,-1,-1,-1,-1,  0, 1, 5, 6, 2,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 5, 6, 3, 3,-1, -1,-1,-1,-1,-1,-1,  0, 1, 5, 6, 3, 3,  5, 6, 3, 0, 1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       4, 6, 2, 1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 4, 6, 2, 1,-1, -1,-1,-1,-1,-1,-1, &
       1, 1, 4, 6, 2,-1, -1,-1,-1,-1,-1,-1,  0, 4, 6, 2,-1,-1,  2, 3, 0, 4, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 4, 6, 3, 3,-1,  4, 6, 3, 0,-1,-1, &
       6, 2, 1, 0, 7,-1, -1,-1,-1,-1,-1,-1,  0, 0, 7, 6, 2, 1,  1, 3, 7, 6, 2,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  0, 0, 7, 6, 2,-1,  2, 3, 7, 6,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  6, 3, 7,-1,-1,-1, &
       6, 7, 3,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 7, 3,-1,-1,  0, 2, 6, 7, 0,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 2, 6, 7, 3,-1,  0, 1, 2, 6, 7, 0, -1,-1,-1,-1,-1,-1,  6, 7, 0, 1, 2,-1, &
       4, 0, 3, 6,-1,-1,  0, 3, 3, 6, 4,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 4, 0, 3,-1,  0, 2, 6, 4,-1,-1, -1,-1,-1,-1,-1,-1,  1, 2, 6, 4, 1,-1, &
      -1,-1,-1,-1,-1,-1,  0, 1, 2, 6, 4,-1, -1,-1,-1,-1,-1,-1,  4, 1, 2, 6,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       5, 1, 0, 3, 6,-1,  0, 3, 3, 6, 5, 1, -1,-1,-1,-1,-1,-1,  1, 3, 3, 6, 5,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1,  0, 2, 6, 5, 1,-1, -1,-1,-1,-1,-1,-1,  1, 2, 6, 5,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  5, 2, 6,-1,-1,-1, &
       5, 7, 3, 2,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 5, 7, 3, 2,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 2, 5, 7, 3,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 5, 7, 3,-1,-1,  0, 1, 5, 7, 0,-1,  0, 1, 5, 7, 3,-1,  5, 7, 0, 1,-1,-1, &
       4, 0, 3, 2, 5,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 2, 5, 4, 0, 3,  0, 2, 2, 5, 4,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 5, 4, 0, 3,-1,  0, 1, 5, 4,-1,-1, -1,-1,-1,-1,-1,-1,  4, 1, 5,-1,-1,-1, &
       4, 7, 3, 2, 1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 1, 4, 7, 3, 2, -1,-1,-1,-1,-1,-1,  0, 4, 7, 3, 2,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       1, 1, 4, 7, 3,-1, -1,-1,-1,-1,-1,-1,  0, 4, 7, 3,-1,-1,  4, 7, 0,-1,-1,-1, &
       0, 3, 2, 1,-1,-1,  0, 3, 3, 2, 1,-1,  0, 0, 3, 2, 1,-1,  1, 3, 3, 2,-1,-1, &
       1, 1, 0, 3, 2,-1, -1,-1,-1,-1,-1,-1,  0, 0, 3, 2,-1,-1,  2, 3, 3,-1,-1,-1, &
       2, 2, 1, 0, 3,-1,  0, 2, 2, 1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 2, 2,-1,-1,-1, &
       1, 1, 0, 3,-1,-1,  0, 1, 1,-1,-1,-1,  0, 0, 3,-1,-1,-1, -1,-1,-1,-1,-1,-1], &
      shape=[6,256])

   integer, dimension(6,0:255), parameter :: hex_poly_v2 = reshape([ &
      -1,-1,-1,-1,-1,-1,  1, 0, 4,-1,-1,-1,  1, 5, 2,-1,-1,-1,  2, 0, 4, 5,-1,-1, &
       2, 6, 3,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 6, 3,-1,-1,  3, 0, 4, 5, 6,-1, &
       3, 7, 0,-1,-1,-1,  1, 3, 7, 4,-1,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 4, 5,-1, &
       2, 6, 7, 0,-1,-1,  1, 2, 6, 7, 4,-1,  1, 5, 6, 7, 0,-1,  4, 5, 6, 7,-1,-1, &
       5, 4, 4,-1,-1,-1,  1, 0, 4, 5,-1,-1, -1,-1,-1,-1,-1,-1,  2, 0, 4, 5, 5,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1,  1, 3, 7, 4, 5,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 4, 5, 5, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  5, 5, 6, 7, 4,-1, &
       5, 6, 5,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 6, 2,-1,-1,  2, 0, 4, 5, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 6, 6, 3,-1,  3, 0, 4, 5, 6, 6, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  5, 6, 6, 7, 4,-1, &
       6, 5, 4, 4,-1,-1,  1, 0, 4, 6, 5,-1,  1, 4, 4, 6, 2,-1,  2, 0, 4, 6,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  3, 0, 4, 6, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 4, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  6, 6, 7, 4,-1,-1, &
       6, 7, 6,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 7, 3,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 6, 7, 3,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 7, 7, 0,-1, -1,-1,-1,-1,-1,-1,  1, 5, 6, 7, 7, 0,  6, 7, 7, 4, 5,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       5, 7, 6, 5,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 7, 6, 2,-1, -1,-1,-1,-1,-1,-1, &
       2, 5, 5, 7, 3,-1, -1,-1,-1,-1,-1,-1,  1, 5, 7, 3,-1,-1,  3, 0, 4, 5, 7,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 5, 7, 7, 0,-1,  5, 7, 7, 4,-1,-1, &
       7, 6, 5, 4, 4,-1, -1,-1,-1,-1,-1,-1,  1, 4, 4, 7, 6, 2,  2, 0, 4, 7, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  1, 4, 4, 7, 3,-1,  3, 0, 4, 7,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  7, 7, 4,-1,-1,-1, &
       7, 4, 7,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       3, 7, 4, 0,-1,-1,  1, 3, 7, 4, 4,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 7, 4, 0,-1,  1, 2, 6, 7, 4, 4, -1,-1,-1,-1,-1,-1,  7, 4, 4, 5, 6,-1, &
       5, 4, 7, 7,-1,-1,  1, 0, 7, 7, 5,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       3, 7, 5, 4, 0,-1,  1, 3, 7, 5,-1,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 5, 5,-1, &
      -1,-1,-1,-1,-1,-1,  1, 2, 6, 7, 5,-1, -1,-1,-1,-1,-1,-1,  5, 5, 6, 7,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       6, 5, 4, 7, 7,-1,  1, 0, 7, 7, 6, 5, -1,-1,-1,-1,-1,-1,  2, 0, 7, 7, 6,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1,  1, 3, 7, 6, 5,-1, -1,-1,-1,-1,-1,-1,  2, 3, 7, 6,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,  6, 6, 7,-1,-1,-1, &
       6, 4, 7, 6,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 4, 7, 3,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       3, 6, 6, 4, 0,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 4, 0,-1,-1,  1, 2, 6, 4, 4,-1,  1, 5, 6, 4, 0,-1,  6, 4, 4, 5,-1,-1, &
       5, 4, 7, 6, 6,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       3, 6, 6, 5, 4, 0,  1, 3, 6, 6, 5,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 6, 5, 4, 0,-1,  1, 2, 6, 5,-1,-1, -1,-1,-1,-1,-1,-1,  5, 5, 6,-1,-1,-1, &
       5, 4, 7, 6, 5,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 5, 5, 4, 7, 3, -1,-1,-1,-1,-1,-1,  1, 5, 4, 7, 3,-1, -1,-1,-1,-1,-1,-1, &
      -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1, &
       2, 5, 5, 4, 0,-1, -1,-1,-1,-1,-1,-1,  1, 5, 4, 0,-1,-1,  5, 4, 4,-1,-1,-1, &
       4, 7, 6, 5,-1,-1,  1, 0, 7, 6, 5,-1,  1, 4, 7, 6, 2,-1,  2, 0, 7, 6,-1,-1, &
       2, 5, 4, 7, 3,-1, -1,-1,-1,-1,-1,-1,  1, 4, 7, 3,-1,-1,  3, 0, 7,-1,-1,-1, &
       3, 6, 5, 4, 0,-1,  1, 3, 6, 5,-1,-1, -1,-1,-1,-1,-1,-1,  2, 3, 6,-1,-1,-1, &
       2, 5, 4, 0,-1,-1,  1, 2, 5,-1,-1,-1,  1, 4, 0,-1,-1,-1, -1,-1,-1,-1,-1,-1], &
      shape=[6,256])

contains

   !> Compute plane distance d such that plane (n,d) cuts cell to give target VF
   !> Cell is axis-aligned box defined by corners lo and hi
   !> Normal n must be unit length
   !> Returns d in global coordinates: n.x = d
   !> Reference: Scardovelli & Zaleski, JCP 164, 228-247 (2000)
   pure function get_plane_dist(normal, lo, hi, VF) result(d)
      implicit none
      real(WP), dimension(3), intent(in) :: normal  !< Unit normal
      real(WP), dimension(3), intent(in) :: lo      !< Lower cell corner
      real(WP), dimension(3), intent(in) :: hi      !< Upper cell corner
      real(WP), intent(in) :: VF                    !< Target volume fraction
      real(WP) :: d
      real(WP) :: m1, m2, m3, mm1, mm2, mm3, tmp
      real(WP) :: norm, factor, alpha, VOFo
      real(WP) :: m12, V1, V2, V3
      real(WP) :: a0, a1, a2, p0, q0, theta
      real(WP), parameter :: eps = 1.0e-15_WP
      real(WP), dimension(3) :: cellsize, cellctr
      
      ! Compute cell size and centroid from corners
      cellsize = hi - lo
      cellctr  = 0.5_WP * (lo + hi)
      
      ! Form scaled normal components
      mm1 = normal(1) * cellsize(1)
      mm2 = normal(2) * cellsize(2)
      mm3 = normal(3) * cellsize(3)
      
      ! Normalize to sum of absolute values = 1
      norm = abs(mm1) + abs(mm2) + abs(mm3)
      if (norm.lt.eps) then
         d = sign(1.0e10_WP, VF - 0.5_WP)
         return
      end if
      mm1 = mm1 / norm
      mm2 = mm2 / norm
      mm3 = mm3 / norm

      
      ! Track offset for negative components
      factor = 0.0_WP
      if (mm1.lt.0.0_WP) factor = factor + mm1
      if (mm2.lt.0.0_WP) factor = factor + mm2
      if (mm3.lt.0.0_WP) factor = factor + mm3
      
      ! Handle VF > 0.5 by symmetry
      VOFo = VF
      if (VF.gt.0.5_WP) VOFo = 1.0_WP - VF
      
      ! === Core S&Z algorithm ===
      ! Take absolute value and sort so that m1 <= m2 <= m3
      m1 = abs(mm1); m2 = abs(mm2); m3 = abs(mm3)
      if (m2.lt.m1) then; tmp = m2; m2 = m1; m1 = tmp; end if
      if (m3.lt.m2) then; tmp = m3; m3 = m2; m2 = tmp; end if
      if (m2.lt.m1) then; tmp = m2; m2 = m1; m1 = tmp; end if
      
      ! Form volume thresholds
      m12 = m1 + m2
      V1 = m1**2 / max(6.0_WP * m2 * m3, eps)
      V2 = V1 + 0.5_WP * (m2 - m1) / max(m3, eps)
      if (m12.le.m3) then
         V3 = 0.5_WP * m12 / max(m3, eps)
      else
         V3 = (m3**2 * (3.0_WP*m12 - m3) + m1**2 * (m1 - 3.0_WP*m3) + &
              m2**2 * (m2 - 3.0_WP*m3)) / max(6.0_WP * m1 * m2 * m3, eps)
      end if
      
      ! Calculate alpha based on regime
      if (VOFo.ge.0.0_WP.and.VOFo.lt.V1) then
         alpha = (6.0_WP * m1 * m2 * m3 * VOFo)**(1.0_WP/3.0_WP)
      else if (VOFo.lt.V2) then
         alpha = 0.5_WP * (m1 + sqrt(m1**2 + 8.0_WP * m2 * m3 * (VOFo - V1)))
      else if (VOFo.lt.V3) then
         a0 = -(m1**3 + m2**3 - 6.0_WP * m1 * m2 * m3 * VOFo)
         a1 = 3.0_WP * (m1**2 + m2**2)
         a2 = -3.0_WP * m12
         p0 = -(a1 / 3.0_WP - a2**2 / 9.0_WP)
         q0 = (a1 * a2 - 3.0_WP * a0) / 6.0_WP - a2**3 / 27.0_WP
         theta = acos(max(-1.0_WP, min(1.0_WP, q0 / sqrt(p0**3)))) / 3.0_WP
         alpha = sqrt(p0) * (sqrt(3.0_WP) * sin(theta) - cos(theta)) - a2 / 3.0_WP
      else
         if (m12.le.m3) then
            alpha = m3 * VOFo + 0.5_WP * m12
         else
            a0 = -0.5_WP * (m1**3 + m2**3 + m3**3 - 6.0_WP * m1 * m2 * m3 * VOFo)
            a1 = 1.5_WP * (m1**2 + m2**2 + m3**2)
            a2 = -1.5_WP
            p0 = -(a1 / 3.0_WP - a2**2 / 9.0_WP)
            q0 = (a1 * a2 - 3.0_WP * a0) / 6.0_WP - a2**3 / 27.0_WP
            theta = acos(max(-1.0_WP, min(1.0_WP, q0 / sqrt(p0**3)))) / 3.0_WP
            alpha = sqrt(p0) * (sqrt(3.0_WP) * sin(theta) - cos(theta)) - a2 / 3.0_WP
         end if
      end if
      ! === End core S&Z algorithm ===
      
      ! Undo VF > 0.5 symmetry
      if (VF.gt.0.5_WP) alpha = 1.0_WP - alpha
      
      ! Adjust for negative normal components
      alpha = alpha + factor
      
      ! Transform from cell-corner origin to cell-centroid origin
      alpha = alpha - 0.5_WP * (mm1 + mm2 + mm3)
      
      ! Scale back to physical coordinates and add centroid offset
      d = alpha * norm + dot_product(normal, cellctr)
      
   end function get_plane_dist
   
   !> Compute volume of tetrahedron given 4 vertices
   !> v(:,1:4) are the vertex coordinates
   pure function tet_vol(v) result(vol)
      implicit none
      real(WP), dimension(3,4), intent(in) :: v
      real(WP) :: vol
      real(WP), dimension(3) :: a,b,c
      a=v(:,1)-v(:,4)
      b=v(:,2)-v(:,4)
      c=v(:,3)-v(:,4)
      vol=(-a(1)*(b(2)*c(3)-c(2)*b(3))&
      &    +a(2)*(b(1)*c(3)-c(1)*b(3))&
      &    -a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
   end function tet_vol

   !> Function that calculates the sign of a tet
   function tet_sign(vert) result(s)
      implicit none
      real(WP) :: s
      real(WP), dimension(3,4), intent(in) :: vert
      real(WP), dimension(3) :: a,b,c
      a=vert(:,1)-vert(:,4); b=vert(:,2)-vert(:,4); c=vert(:,3)-vert(:,4)
      s=sign(1.0_WP,-(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2))))
   end function tet_sign
   
   !> Cut a tetrahedron by a plane and return liquid/gas volumes and barycenters
   !> Input:  v(:,1:4) = 4 tet vertices, plane(1:4) = [nx,ny,nz,d] where n.x=d defines plane
   !> Output: vol_liq, vol_gas = phase volumes
   !>         bary_liq(3), bary_gas(3) = volume-weighted barycenters
   pure subroutine cut_tet_vol(v, plane, vol_liq, vol_gas, bary_liq, bary_gas)
      implicit none
      real(WP), dimension(3,4), intent(in)  :: v
      real(WP), dimension(4),   intent(in)  :: plane
      real(WP),                 intent(out) :: vol_liq, vol_gas
      real(WP), dimension(3),   intent(out) :: bary_liq, bary_gas
      real(WP), dimension(4)   :: d
      real(WP), dimension(3,8) :: vert
      real(WP), dimension(3)   :: a, b, c, bary
      real(WP) :: denom, mu, my_vol
      integer  :: icase, n1, v1, v2
      
      ! Compute signed distance of each vertex to plane
      d(1) = plane(1)*v(1,1) + plane(2)*v(2,1) + plane(3)*v(3,1) - plane(4)
      d(2) = plane(1)*v(1,2) + plane(2)*v(2,2) + plane(3)*v(3,2) - plane(4)
      d(3) = plane(1)*v(1,3) + plane(2)*v(2,3) + plane(3)*v(3,3) - plane(4)
      d(4) = plane(1)*v(1,4) + plane(2)*v(2,4) + plane(3)*v(3,4) - plane(4)
      
      ! Determine cutting case (1-16 based on which vertices are above plane)
      icase = 1 + int(0.5_WP+sign(0.5_WP,d(1))) + 2*int(0.5_WP+sign(0.5_WP,d(2))) &
                + 4*int(0.5_WP+sign(0.5_WP,d(3))) + 8*int(0.5_WP+sign(0.5_WP,d(4)))
      
      ! Copy original vertices
      vert(:,1:4) = v(:,1:4)
      
      ! Create interpolated vertices on cut plane
      do n1 = 1, cut_nvert(icase)
         v1 = cut_v1(n1,icase); v2 = cut_v2(n1,icase)
         denom = d(v2) - d(v1)
         if (abs(denom).ge.tiny(1.0_WP)) then
            mu = max(0.0_WP, min(1.0_WP, -d(v1)/denom))
         else
            mu = 0.0_WP
         end if
         vert(:,4+n1) = (1.0_WP-mu)*vert(:,v1) + mu*vert(:,v2)
      end do
      
      ! Initialize outputs
      vol_liq = 0.0_WP; vol_gas = 0.0_WP
      bary_liq = 0.0_WP; bary_gas = 0.0_WP
      
      ! Gas tets: from 1 to cut_nntet-1
      do n1 = 1, cut_nntet(icase)-1
         a = vert(:,cut_vtet(1,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         b = vert(:,cut_vtet(2,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         c = vert(:,cut_vtet(3,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         my_vol = abs(a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
         bary = 0.25_WP * (vert(:,cut_vtet(1,n1,icase)) + vert(:,cut_vtet(2,n1,icase)) &
                         + vert(:,cut_vtet(3,n1,icase)) + vert(:,cut_vtet(4,n1,icase)))
         vol_gas = vol_gas + my_vol
         bary_gas = bary_gas + my_vol * bary
      end do
      
      ! Liquid tets: from cut_ntets down to cut_nntet
      do n1 = cut_ntets(icase), cut_nntet(icase), -1
         a = vert(:,cut_vtet(1,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         b = vert(:,cut_vtet(2,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         c = vert(:,cut_vtet(3,n1,icase)) - vert(:,cut_vtet(4,n1,icase))
         my_vol = abs(a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
         bary = 0.25_WP * (vert(:,cut_vtet(1,n1,icase)) + vert(:,cut_vtet(2,n1,icase)) &
                         + vert(:,cut_vtet(3,n1,icase)) + vert(:,cut_vtet(4,n1,icase)))
         vol_liq = vol_liq + my_vol
         bary_liq = bary_liq + my_vol * bary
      end do
      
      ! Normalize barycenters
      if (vol_liq.gt.tiny(1.0_WP)) bary_liq = bary_liq / vol_liq
      if (vol_gas.gt.tiny(1.0_WP)) bary_gas = bary_gas / vol_gas
      
   end subroutine cut_tet_vol
   
   !> Compute signed volume of flux polyhedron
   !> face(:,1:8) = 8 vertices (4 at time t, 4 back-projected)
   !> Decomposes into 6 tets using tet_map and sums signed volumes
   pure function flux_polyhedron_vol(face) result(vol)
      implicit none
      real(WP), dimension(3,8), intent(in) :: face
      real(WP) :: vol
      real(WP), dimension(3) :: a, b, c
      integer :: ntet
      vol = 0.0_WP
      do ntet = 1, 6
         a = face(:,tet_map(1,ntet)) - face(:,tet_map(4,ntet))
         b = face(:,tet_map(2,ntet)) - face(:,tet_map(4,ntet))
         c = face(:,tet_map(3,ntet)) - face(:,tet_map(4,ntet))
         vol = vol + (a(1)*(b(2)*c(3)-c(2)*b(3)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2))) / 6.0_WP
      end do
   end function flux_polyhedron_vol

   !> Adjust flux polyhedron to enforce target volume
   !> Uses IRL's direction-independent approach: moves vertex 9 along back-face normal
   subroutine correct_flux_poly(poly,target_volume)
      implicit none
      real(WP), dimension(3,9), intent(inout) :: poly
      real(WP), intent(in) :: target_volume
      real(WP) :: starting_volume,needed_change,mag,adjustment
      real(WP), dimension(3) :: cross_sum,dir
      real(WP), dimension(3) :: e1,e2,c1,c2
      integer :: n
      ! Compute starting volume
      starting_volume=0.0_WP
      do n=1,8; starting_volume=starting_volume+tet_vol([poly(:,tet_map(1,n)),poly(:,tet_map(2,n)),poly(:,tet_map(3,n)),poly(:,tet_map(4,n))]); end do
      needed_change=target_volume-starting_volume
      ! Compute volume gradient for tets 7 and 8 (the only tets using vertex 9)
      ! For tet (A,B,C,D), gradient w.r.t. D is -(1/6)*((B-A)×(C-A))
      ! Tet 7 = (5,6,8,9): gradient = -((v6-v5)×(v8-v5)). Take edges from vertex 5
      e1=poly(:,6)-poly(:,5); e2=poly(:,8)-poly(:,5)
      c1=[e1(2)*e2(3)-e1(3)*e2(2),e1(3)*e2(1)-e1(1)*e2(3),e1(1)*e2(2)-e1(2)*e2(1)]
      ! Tet 8 = (6,7,8,9): gradient = -((v7-v6)×(v8-v6)). Take edges from vertex 6
      e1=poly(:,7)-poly(:,6); e2=poly(:,8)-poly(:,6)
      c2=[e1(2)*e2(3)-e1(3)*e2(2),e1(3)*e2(1)-e1(1)*e2(3),e1(1)*e2(2)-e1(2)*e2(1)]
      ! Total gradient (negative sign absorbed into adjustment formula)
      cross_sum=c1+c2
      ! Compute adjustment along normal direction
      mag=sqrt(cross_sum(1)**2+cross_sum(2)**2+cross_sum(3)**2)
      adjustment=6.0_WP*needed_change/max(mag,tiny(1.0_WP))
      dir=cross_sum/max(mag,tiny(1.0_WP))
      ! Move vertex 9
      poly(:,9)=poly(:,9)+adjustment*dir
   end subroutine correct_flux_poly

   !> Compute total signed volume and volume-weighted barycenter of a
   !> corrected 9-vertex flux polyhedron using tet_map(4,8) decomposition
   !> No PLIC cutting, no grid-plane recursion.
   pure subroutine flux_poly_moments(poly,vol,bary)
      implicit none
      real(WP), dimension(3,9), intent(in)  :: poly
      real(WP),                 intent(out) :: vol
      real(WP), dimension(3),   intent(out) :: bary
      real(WP), dimension(3) :: a,b,c,centroid
      real(WP) :: svol
      integer :: n
      vol =0.0_WP
      bary=0.0_WP
      do n=1,8
         a=poly(:,tet_map(1,n))-poly(:,tet_map(4,n))
         b=poly(:,tet_map(2,n))-poly(:,tet_map(4,n))
         c=poly(:,tet_map(3,n))-poly(:,tet_map(4,n))
         svol=(-a(1)*(b(2)*c(3)-c(2)*b(3))+a(2)*(b(1)*c(3)-c(1)*b(3))-a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
         centroid=0.25_WP*(poly(:,tet_map(1,n))+poly(:,tet_map(2,n))+poly(:,tet_map(3,n))+poly(:,tet_map(4,n)))
         vol=vol+svol
         bary=bary+svol*centroid
      end do
   end subroutine flux_poly_moments

   !> Cut a hex cell by a plane and compute liquid/gas volumes and barycenters
   !> hex(:,1:8) = 8 vertices of hex cell (standard ordering)
   !> plane(1:3) = normal, plane(4) = distance (n·x = d)
   !> Returns vol_liq/vol_gas and bary_liq/bary_gas
   subroutine cut_hex_vol(hex, plane, vol_liq, vol_gas, bary_liq, bary_gas)
      implicit none
      real(WP), dimension(3,8), intent(in) :: hex
      real(WP), dimension(4), intent(in) :: plane
      real(WP), intent(out) :: vol_liq, vol_gas
      real(WP), dimension(3), intent(out) :: bary_liq, bary_gas
      ! Local
      real(WP), dimension(3,4) :: tet
      real(WP), dimension(3) :: center
      real(WP) :: vl, vg
      real(WP), dimension(3) :: bl, bg
      integer :: ntet, v
      ! Hex->tet decomposition: 6 tets using center point
      ! hex_tet_map: indices into hex vertices for 6 tets (each shares center)
      integer, dimension(3,6), parameter :: hex_tet_map = reshape([ &
         1, 2, 3,  1, 3, 4,  1, 2, 6,  1, 6, 5, &
         5, 6, 7,  5, 7, 8], shape(hex_tet_map))
      integer, dimension(3,6), parameter :: hex_tet_top = reshape([ &
         5, 6, 7,  5, 7, 8,  4, 3, 7,  4, 7, 8, &
         1, 2, 3,  1, 3, 4], shape(hex_tet_top))
      ! Actually use standard 5-tet decomposition of hex
      integer, dimension(4,5), parameter :: tet5 = reshape([ &
         1, 2, 4, 5,  2, 3, 4, 7,  2, 5, 6, 7, &
         4, 5, 7, 8,  2, 4, 5, 7], shape(tet5))
      
      vol_liq = 0.0_WP; vol_gas = 0.0_WP
      bary_liq = 0.0_WP; bary_gas = 0.0_WP
      
      ! Decompose hex into 5 tetrahedra
      do ntet = 1, 5
         ! Build tetrahedron
         do v = 1, 4
            tet(:,v) = hex(:,tet5(v,ntet))
         end do
         ! Cut tet by plane
         call cut_tet_vol(tet, plane, vl, vg, bl, bg)
         ! Accumulate volumes and volume-weighted barycenters
         bary_liq = bary_liq + vl * bl
         bary_gas = bary_gas + vg * bg
         vol_liq = vol_liq + vl
         vol_gas = vol_gas + vg
      end do
      
      ! Normalize barycenters
      if (vol_liq.gt.tiny(1.0_WP)) bary_liq = bary_liq / vol_liq
      if (vol_gas.gt.tiny(1.0_WP)) bary_gas = bary_gas / vol_gas
      
   end subroutine cut_hex_vol

   !> =========================================================================
   !> Extract PLIC polygon from hex cell using IRL's algorithm
   !> Given 8 hex vertices and a plane (normal, distance), compute the 
   !> intersection polygon vertices.
   !> 
   !> Hex vertex ordering (IRL convention):
   !>   0:(+,-,-), 1:(+,+,-), 2:(+,+,+), 3:(+,-,+)
   !>   4:(-,-,-), 5:(-,+,-), 6:(-,+,+), 7:(-,-,+)
   !>
   !> For axis-aligned cell [xlo,ylo,zlo] to [xhi,yhi,zhi]:
   !>   hex(:,1) = [xhi, ylo, zlo]  ! vertex 0
   !>   hex(:,2) = [xhi, yhi, zlo]  ! vertex 1
   !>   hex(:,3) = [xhi, yhi, zhi]  ! vertex 2
   !>   hex(:,4) = [xhi, ylo, zhi]  ! vertex 3
   !>   hex(:,5) = [xlo, ylo, zlo]  ! vertex 4
   !>   hex(:,6) = [xlo, yhi, zlo]  ! vertex 5
   !>   hex(:,7) = [xlo, yhi, zhi]  ! vertex 6
   !>   hex(:,8) = [xlo, ylo, zhi]  ! vertex 7
   !>
   !> Returns:
   !>   nvert: number of polygon vertices (0-6)
   !>   poly(:,1:nvert): polygon vertices in cyclic order
   !> =========================================================================
   pure subroutine cut_hex_polygon(hex, plane, nvert, poly)
      implicit none
      real(WP), dimension(3,8), intent(in)  :: hex      !< 8 hex vertices
      real(WP), dimension(4),   intent(in)  :: plane    !< (nx, ny, nz, d)
      integer,                  intent(out) :: nvert    !< Number of polygon vertices
      real(WP), dimension(3,6), intent(out) :: poly     !< Polygon vertices (up to 6)
      
      ! Local variables
      real(WP), dimension(8) :: dist       ! Signed distance to plane for each vertex
      integer :: icase                     ! Cutting case (0-255)
      integer :: n, v1, v2
      real(WP) :: mu, denom
      
      ! Initialize output
      nvert = 0
      poly = 0.0_WP
      
      ! Compute signed distance from plane for each vertex
      ! Plane equation: n.x - d = 0, so distance = n.x - d
      do n = 1, 8
         dist(n) = plane(1)*hex(1,n) + plane(2)*hex(2,n) + plane(3)*hex(3,n) - plane(4)
      end do
      
      ! Determine cutting case from sign pattern (0-indexed vertices)
      ! Case bit i is set if vertex i is above plane (dist > 0)
      icase = 0
      if (dist(1).gt.0.0_WP) icase = icase + 1    ! vertex 0
      if (dist(2).gt.0.0_WP) icase = icase + 2    ! vertex 1
      if (dist(3).gt.0.0_WP) icase = icase + 4    ! vertex 2
      if (dist(4).gt.0.0_WP) icase = icase + 8    ! vertex 3
      if (dist(5).gt.0.0_WP) icase = icase + 16   ! vertex 4
      if (dist(6).gt.0.0_WP) icase = icase + 32   ! vertex 5
      if (dist(7).gt.0.0_WP) icase = icase + 64   ! vertex 6
      if (dist(8).gt.0.0_WP) icase = icase + 128  ! vertex 7
      
      ! Get number of polygon vertices for this case
      nvert = hex_poly_nvert(icase)
      
      ! Handle no intersection or ambiguous cases
      if (nvert.le.0) then
         nvert = 0
         return
      end if
      
      ! Compute intersection points along cut edges
      do n = 1, nvert
         ! Get edge endpoints (0-indexed in tables, convert to 1-indexed)
         v1 = hex_poly_v1(n, icase) + 1
         v2 = hex_poly_v2(n, icase) + 1
         
         ! Linear interpolation: find point where dist = 0
         ! mu = -dist(v1) / (dist(v2) - dist(v1))
         denom = dist(v2) - dist(v1)
         if (abs(denom).ge.tiny(1.0_WP)) then
            mu = max(0.0_WP, min(1.0_WP, -dist(v1)/denom))
         else
            mu = 0.5_WP
         end if
         
         ! Interpolated point
         poly(:,n) = (1.0_WP - mu)*hex(:,v1) + mu*hex(:,v2)
      end do
      
   end subroutine cut_hex_polygon
   
   !> =========================================================================
   !> Cheap vertex count for hex-plane intersection (no interpolation)
   !> Returns the number of polygon vertices from cutting case lookup only.
   !> This is much cheaper than cut_hex_polygon since it skips interpolation.
   !> =========================================================================
   pure function get_hex_poly_nvert(hex, plane) result(nvert)
      implicit none
      real(WP), dimension(3,8), intent(in) :: hex      !< 8 hex vertices
      real(WP), dimension(4),   intent(in) :: plane    !< (nx, ny, nz, d)
      integer :: nvert
      
      real(WP), dimension(8) :: dist
      integer :: icase, n
      
      ! Compute signed distance from plane for each vertex
      do n = 1, 8
         dist(n) = plane(1)*hex(1,n) + plane(2)*hex(2,n) + plane(3)*hex(3,n) - plane(4)
      end do
      
      ! Determine cutting case from sign pattern
      icase = 0
      if (dist(1).gt.0.0_WP) icase = icase + 1
      if (dist(2).gt.0.0_WP) icase = icase + 2
      if (dist(3).gt.0.0_WP) icase = icase + 4
      if (dist(4).gt.0.0_WP) icase = icase + 8
      if (dist(5).gt.0.0_WP) icase = icase + 16
      if (dist(6).gt.0.0_WP) icase = icase + 32
      if (dist(7).gt.0.0_WP) icase = icase + 64
      if (dist(8).gt.0.0_WP) icase = icase + 128
      
      ! Get vertex count from lookup table
      nvert = hex_poly_nvert(icase)
      if (nvert.lt.0) nvert = 0
      
   end function get_hex_poly_nvert

   !> =========================================================================
   !> Remap all nodes in a box via RK2 back-projection (staggered velocity)
   !>
   !> For each node (ii,jj,kk) in [bx%lo:bx%hi+1]:
   !>   Stage 1: velocity at vertex
   !>   Stage 2: velocity at half-step position
   !>
   !> proj(ii,jj,kk,1:3) = projected position of node (ii,jj,kk)
   !> =========================================================================
   pure subroutine remap_box_staggered(bx,dt,xlo,ylo,zlo,dx,dy,dz,U,Ulo,V,Vlo,W,Wlo,proj,Projlo)
      use amrex_amr_module, only: amrex_box
      implicit none
      type(amrex_box), intent(in) :: bx
      real(WP), intent(in) :: dt
      real(WP), intent(in) :: xlo,ylo,zlo
      real(WP), intent(in) :: dx,dy,dz
      integer, dimension(4), intent(in) :: Ulo,Vlo,Wlo
      real(WP), dimension(Ulo(1):,Ulo(2):,Ulo(3):,Ulo(4):), intent(in) :: U
      real(WP), dimension(Vlo(1):,Vlo(2):,Vlo(3):,Vlo(4):), intent(in) :: V
      real(WP), dimension(Wlo(1):,Wlo(2):,Wlo(3):,Wlo(4):), intent(in) :: W
      integer, dimension(4), intent(in) :: Projlo
      real(WP), dimension(Projlo(1):,Projlo(2):,Projlo(3):,Projlo(4):), intent(out) :: proj
      ! Local
      integer :: ii,jj,kk
      real(WP) :: px,py,pz
      real(WP) :: u1,v1,w1
      real(WP) :: hx,hy,hz
      real(WP) :: u2,v2,w2
      real(WP) :: dxi,dyi,dzi
      ! For stage 2 trilinear
      integer  :: ipc,jpc,kpc
      integer  :: ipu,jpv,kpw
      real(WP) :: wxc,wyc,wzc,wxc2,wyc2,wzc2
      real(WP) :: wxu,wyv,wzw,wxu2,wyv2,wzw2
      ! Precompute inverse mesh
      dxi=1.0_WP/dx; dyi=1.0_WP/dy; dzi=1.0_WP/dz
      ! Loop over all nodes in box
      do kk=bx%lo(3),bx%hi(3)
         do jj=bx%lo(2),bx%hi(2)
            do ii=bx%lo(1),bx%hi(1)
               ! Physical position of node (ii,jj,kk)
               px=xlo+real(ii,WP)*dx
               py=ylo+real(jj,WP)*dy
               pz=zlo+real(kk,WP)*dz
               ! Stage 1
               hx=px-0.5_WP*dt*0.25_WP*sum(U(ii,jj-1:jj,kk-1:kk,1))
               hy=py-0.5_WP*dt*0.25_WP*sum(V(ii-1:ii,jj,kk-1:kk,1))
               hz=pz-0.5_WP*dt*0.25_WP*sum(W(ii-1:ii,jj-1:jj,kk,1))
               ! Stage 2
               ! Cell-centered indices
               ipc = floor((hx - xlo)*dxi - 0.5_WP)
               jpc = floor((hy - ylo)*dyi - 0.5_WP)
               kpc = floor((hz - zlo)*dzi - 0.5_WP)
               ! Face-centered indices
               ipu = floor((hx - xlo)*dxi)
               jpv = floor((hy - ylo)*dyi)
               kpw = floor((hz - zlo)*dzi)
               ! Clamp to array bounds
               ipu = max(lbound(U,1), min(ubound(U,1)-1, ipu))
               jpc = max(lbound(U,2), min(ubound(U,2)-1, jpc))
               kpc = max(lbound(U,3), min(ubound(U,3)-1, kpc))
               ipc = max(lbound(V,1), min(ubound(V,1)-1, ipc))
               jpv = max(lbound(V,2), min(ubound(V,2)-1, jpv))
               kpw = max(lbound(W,3), min(ubound(W,3)-1, kpw))
               ! Cell-centered weights
               wxc = (hx - (xlo + (real(ipc,WP)+0.5_WP)*dx))*dxi
               wyc = (hy - (ylo + (real(jpc,WP)+0.5_WP)*dy))*dyi
               wzc = (hz - (zlo + (real(kpc,WP)+0.5_WP)*dz))*dzi
               wxc = max(0.0_WP, min(1.0_WP, wxc)); wxc2 = 1.0_WP - wxc
               wyc = max(0.0_WP, min(1.0_WP, wyc)); wyc2 = 1.0_WP - wyc
               wzc = max(0.0_WP, min(1.0_WP, wzc)); wzc2 = 1.0_WP - wzc
               ! Face-centered weights
               wxu = (hx - (xlo + real(ipu,WP)*dx))*dxi
               wyv = (hy - (ylo + real(jpv,WP)*dy))*dyi
               wzw = (hz - (zlo + real(kpw,WP)*dz))*dzi
               wxu = max(0.0_WP, min(1.0_WP, wxu)); wxu2 = 1.0_WP - wxu
               wyv = max(0.0_WP, min(1.0_WP, wyv)); wyv2 = 1.0_WP - wyv
               wzw = max(0.0_WP, min(1.0_WP, wzw)); wzw2 = 1.0_WP - wzw
               ! U at x-faces: face-centered in x, cell-centered in y,z
               u2 = wzc *(wyc *(wxu *U(ipu+1,jpc+1,kpc+1,1)+wxu2*U(ipu,jpc+1,kpc+1,1)) + &
               &          wyc2*(wxu *U(ipu+1,jpc  ,kpc+1,1)+wxu2*U(ipu,jpc  ,kpc+1,1))) + &
               &    wzc2*(wyc *(wxu *U(ipu+1,jpc+1,kpc  ,1)+wxu2*U(ipu,jpc+1,kpc  ,1)) + &
               &          wyc2*(wxu *U(ipu+1,jpc  ,kpc  ,1)+wxu2*U(ipu,jpc  ,kpc  ,1)))
               ! V at y-faces: cell-centered in x, face-centered in y, cell-centered in z
               v2 = wzc *(wyv *(wxc *V(ipc+1,jpv+1,kpc+1,1)+wxc2*V(ipc,jpv+1,kpc+1,1)) + &
               &          wyv2*(wxc *V(ipc+1,jpv  ,kpc+1,1)+wxc2*V(ipc,jpv  ,kpc+1,1))) + &
               &    wzc2*(wyv *(wxc *V(ipc+1,jpv+1,kpc  ,1)+wxc2*V(ipc,jpv+1,kpc  ,1)) + &
               &          wyv2*(wxc *V(ipc+1,jpv  ,kpc  ,1)+wxc2*V(ipc,jpv  ,kpc  ,1)))
               ! W at z-faces: cell-centered in x,y, face-centered in z
               w2 = wzw *(wyc *(wxc *W(ipc+1,jpc+1,kpw+1,1)+wxc2*W(ipc,jpc+1,kpw+1,1)) + &
               &          wyc2*(wxc *W(ipc+1,jpc  ,kpw+1,1)+wxc2*W(ipc,jpc  ,kpw+1,1))) + &
               &    wzw2*(wyc *(wxc *W(ipc+1,jpc+1,kpw  ,1)+wxc2*W(ipc,jpc+1,kpw  ,1)) + &
               &          wyc2*(wxc *W(ipc+1,jpc  ,kpw  ,1)+wxc2*W(ipc,jpc  ,kpw  ,1)))
               ! Full-step projected position
               proj(ii,jj,kk,1) = px - dt*u2
               proj(ii,jj,kk,2) = py - dt*v2
               proj(ii,jj,kk,3) = pz - dt*w2
            end do
         end do
      end do
   end subroutine remap_box_staggered

   !> =========================================================================
   !> Remap all nodes in a box via RK2 back-projection (collocated velocity)
   !>
   !> For each node (ii,jj,kk) in [bx%lo:bx%hi+1]:
   !>   Stage 1: velocity at vertex from 8 surrounding cell centers
   !>   Stage 2: velocity at half-step position (cell-centered trilinear)
   !>
   !> proj(ii,jj,kk,1:3) = projected position of node (ii,jj,kk)
   !> =========================================================================
   pure subroutine remap_box_collocated(bx,dt,xlo,ylo,zlo,dx,dy,dz,U,Ulo,V,Vlo,W,Wlo,proj,Projlo)
      use amrex_amr_module, only: amrex_box
      implicit none
      type(amrex_box), intent(in) :: bx
      real(WP), intent(in) :: dt
      real(WP), intent(in) :: xlo,ylo,zlo
      real(WP), intent(in) :: dx,dy,dz
      integer, dimension(4), intent(in) :: Ulo,Vlo,Wlo
      real(WP), dimension(Ulo(1):,Ulo(2):,Ulo(3):,Ulo(4):), intent(in) :: U
      real(WP), dimension(Vlo(1):,Vlo(2):,Vlo(3):,Vlo(4):), intent(in) :: V
      real(WP), dimension(Wlo(1):,Wlo(2):,Wlo(3):,Wlo(4):), intent(in) :: W
      integer, dimension(4), intent(in) :: Projlo
      real(WP), dimension(Projlo(1):,Projlo(2):,Projlo(3):,Projlo(4):), intent(out) :: proj
      ! Local
      integer :: ii,jj,kk
      real(WP) :: px,py,pz
      real(WP) :: hx,hy,hz
      real(WP) :: u2,v2,w2
      real(WP) :: dxi,dyi,dzi
      ! For trilinear interpolation
      integer  :: ipc,jpc,kpc
      real(WP) :: wxc,wyc,wzc,wxc2,wyc2,wzc2
      ! Precompute inverse mesh
      dxi=1.0_WP/dx; dyi=1.0_WP/dy; dzi=1.0_WP/dz
      ! Loop over all nodes in box
      do kk=bx%lo(3),bx%hi(3)
         do jj=bx%lo(2),bx%hi(2)
            do ii=bx%lo(1),bx%hi(1)
               ! Physical position of node (ii,jj,kk)
               px=xlo+real(ii,WP)*dx
               py=ylo+real(jj,WP)*dy
               pz=zlo+real(kk,WP)*dz
               ! Stage 1: average of 8 neighboring cell centers
               hx=px-0.5_WP*dt*0.125_WP*sum(U(ii-1:ii,jj-1:jj,kk-1:kk,1))
               hy=py-0.5_WP*dt*0.125_WP*sum(V(ii-1:ii,jj-1:jj,kk-1:kk,1))
               hz=pz-0.5_WP*dt*0.125_WP*sum(W(ii-1:ii,jj-1:jj,kk-1:kk,1))
               ! Stage 2: cell-centered trilinear interpolation
               ipc = floor((hx - xlo)*dxi - 0.5_WP)
               jpc = floor((hy - ylo)*dyi - 0.5_WP)
               kpc = floor((hz - zlo)*dzi - 0.5_WP)
               ! Clamp to array bounds
               ipc = max(lbound(U,1), min(ubound(U,1)-1, ipc))
               jpc = max(lbound(U,2), min(ubound(U,2)-1, jpc))
               kpc = max(lbound(U,3), min(ubound(U,3)-1, kpc))
               ! Cell-centered weights
               wxc = (hx - (xlo + (real(ipc,WP)+0.5_WP)*dx))*dxi
               wyc = (hy - (ylo + (real(jpc,WP)+0.5_WP)*dy))*dyi
               wzc = (hz - (zlo + (real(kpc,WP)+0.5_WP)*dz))*dzi
               wxc = max(0.0_WP, min(1.0_WP, wxc)); wxc2 = 1.0_WP - wxc
               wyc = max(0.0_WP, min(1.0_WP, wyc)); wyc2 = 1.0_WP - wyc
               wzc = max(0.0_WP, min(1.0_WP, wzc)); wzc2 = 1.0_WP - wzc
               ! U (cell-centered)
               u2 = wzc *(wyc *(wxc *U(ipc+1,jpc+1,kpc+1,1)+wxc2*U(ipc,jpc+1,kpc+1,1)) + &
               &          wyc2*(wxc *U(ipc+1,jpc  ,kpc+1,1)+wxc2*U(ipc,jpc  ,kpc+1,1))) + &
               &    wzc2*(wyc *(wxc *U(ipc+1,jpc+1,kpc  ,1)+wxc2*U(ipc,jpc+1,kpc  ,1)) + &
               &          wyc2*(wxc *U(ipc+1,jpc  ,kpc  ,1)+wxc2*U(ipc,jpc  ,kpc  ,1)))
               ! V (cell-centered)
               v2 = wzc *(wyc *(wxc *V(ipc+1,jpc+1,kpc+1,1)+wxc2*V(ipc,jpc+1,kpc+1,1)) + &
               &          wyc2*(wxc *V(ipc+1,jpc  ,kpc+1,1)+wxc2*V(ipc,jpc  ,kpc+1,1))) + &
               &    wzc2*(wyc *(wxc *V(ipc+1,jpc+1,kpc  ,1)+wxc2*V(ipc,jpc+1,kpc  ,1)) + &
               &          wyc2*(wxc *V(ipc+1,jpc  ,kpc  ,1)+wxc2*V(ipc,jpc  ,kpc  ,1)))
               ! W (cell-centered)
               w2 = wzc *(wyc *(wxc *W(ipc+1,jpc+1,kpc+1,1)+wxc2*W(ipc,jpc+1,kpc+1,1)) + &
               &          wyc2*(wxc *W(ipc+1,jpc  ,kpc+1,1)+wxc2*W(ipc,jpc  ,kpc+1,1))) + &
               &    wzc2*(wyc *(wxc *W(ipc+1,jpc+1,kpc  ,1)+wxc2*W(ipc,jpc+1,kpc  ,1)) + &
               &          wyc2*(wxc *W(ipc+1,jpc  ,kpc  ,1)+wxc2*W(ipc,jpc  ,kpc  ,1)))
               ! Full-step projected position
               proj(ii,jj,kk,1) = px - dt*u2
               proj(ii,jj,kk,2) = py - dt*v2
               proj(ii,jj,kk,3) = pz - dt*w2
            end do
         end do
      end do
   end subroutine remap_box_collocated

   !> RK2 vertex projection using collocated velocity
   pure function project_rk2_collocated(p1,dt,pU,Ulo,pV,Vlo,pW,Wlo,xlo,ylo,zlo,dxi,dyi,dzi) result(p2)
      implicit none
      real(WP), dimension(3), intent(in) :: p1
      real(WP), intent(in) :: dt,xlo,ylo,zlo,dxi,dyi,dzi
      integer, dimension(4), intent(in) :: Ulo,Vlo,Wlo
      real(WP), dimension(Ulo(1):,Ulo(2):,Ulo(3):,Ulo(4):), contiguous, intent(in) :: pU
      real(WP), dimension(Vlo(1):,Vlo(2):,Vlo(3):,Vlo(4):), contiguous, intent(in) :: pV
      real(WP), dimension(Wlo(1):,Wlo(2):,Wlo(3):,Wlo(4):), contiguous, intent(in) :: pW
      real(WP), dimension(3) :: p2,pm,vel
      integer :: ipc,jpc,kpc
      real(WP) :: wx1,wy1,wz1,wx2,wy2,wz2
      integer :: s
      ! Two-stage RK2
      do s=1,2
         ! Pick interpolation point
         if (s.eq.1) then; pm=p1; else; pm=0.5_WP*(p1+p2); end if
         ! Cell-centered indices
         ipc=floor((pm(1)-xlo)*dxi-0.5_WP)
         jpc=floor((pm(2)-ylo)*dyi-0.5_WP)
         kpc=floor((pm(3)-zlo)*dzi-0.5_WP)
         ! Clamp to array bounds
         ipc=max(lbound(pU,1),min(ubound(pU,1)-1,ipc))
         jpc=max(lbound(pU,2),min(ubound(pU,2)-1,jpc))
         kpc=max(lbound(pU,3),min(ubound(pU,3)-1,kpc))
         ! Weights
         wx1=(pm(1)-(xlo+(real(ipc,WP)+0.5_WP)/dxi))*dxi; wx1=max(0.0_WP,min(1.0_WP,wx1)); wx2=1.0_WP-wx1
         wy1=(pm(2)-(ylo+(real(jpc,WP)+0.5_WP)/dyi))*dyi; wy1=max(0.0_WP,min(1.0_WP,wy1)); wy2=1.0_WP-wy1
         wz1=(pm(3)-(zlo+(real(kpc,WP)+0.5_WP)/dzi))*dzi; wz1=max(0.0_WP,min(1.0_WP,wz1)); wz2=1.0_WP-wz1
         ! Trilinear interpolation (all cell-centered)
         vel(1)=wz1*(wy1*(wx1*pU(ipc+1,jpc+1,kpc+1,1)+wx2*pU(ipc,jpc+1,kpc+1,1))+ &
         &           wy2*(wx1*pU(ipc+1,jpc  ,kpc+1,1)+wx2*pU(ipc,jpc  ,kpc+1,1)))+&
         &      wz2*(wy1*(wx1*pU(ipc+1,jpc+1,kpc  ,1)+wx2*pU(ipc,jpc+1,kpc  ,1))+ &
         &           wy2*(wx1*pU(ipc+1,jpc  ,kpc  ,1)+wx2*pU(ipc,jpc  ,kpc  ,1)))
         vel(2)=wz1*(wy1*(wx1*pV(ipc+1,jpc+1,kpc+1,1)+wx2*pV(ipc,jpc+1,kpc+1,1))+ &
         &           wy2*(wx1*pV(ipc+1,jpc  ,kpc+1,1)+wx2*pV(ipc,jpc  ,kpc+1,1)))+&
         &      wz2*(wy1*(wx1*pV(ipc+1,jpc+1,kpc  ,1)+wx2*pV(ipc,jpc+1,kpc  ,1))+ &
         &           wy2*(wx1*pV(ipc+1,jpc  ,kpc  ,1)+wx2*pV(ipc,jpc  ,kpc  ,1)))
         vel(3)=wz1*(wy1*(wx1*pW(ipc+1,jpc+1,kpc+1,1)+wx2*pW(ipc,jpc+1,kpc+1,1))+ &
         &           wy2*(wx1*pW(ipc+1,jpc  ,kpc+1,1)+wx2*pW(ipc,jpc  ,kpc+1,1)))+&
         &      wz2*(wy1*(wx1*pW(ipc+1,jpc+1,kpc  ,1)+wx2*pW(ipc,jpc+1,kpc  ,1))+ &
         &           wy2*(wx1*pW(ipc+1,jpc  ,kpc  ,1)+wx2*pW(ipc,jpc  ,kpc  ,1)))
         ! Project
         p2=p1+dt*vel
      end do
   end function project_rk2_collocated

   !> RK2 vertex projection using staggered velocity
   pure function project_rk2_staggered(p1,dt,pU,Ulo,pV,Vlo,pW,Wlo,xlo,ylo,zlo,dxi,dyi,dzi) result(p2)
      implicit none
      real(WP), dimension(3), intent(in) :: p1
      real(WP), intent(in) :: dt,xlo,ylo,zlo,dxi,dyi,dzi
      integer, dimension(4), intent(in) :: Ulo,Vlo,Wlo
      real(WP), dimension(Ulo(1):,Ulo(2):,Ulo(3):,Ulo(4):), contiguous, intent(in) :: pU
      real(WP), dimension(Vlo(1):,Vlo(2):,Vlo(3):,Vlo(4):), contiguous, intent(in) :: pV
      real(WP), dimension(Wlo(1):,Wlo(2):,Wlo(3):,Wlo(4):), contiguous, intent(in) :: pW
      real(WP), dimension(3) :: p2,pm,vel
      integer :: ipc,jpc,kpc,ipu,jpv,kpw
      real(WP) :: wxc1,wyc1,wzc1,wxc2,wyc2,wzc2
      real(WP) :: wxu1,wyv1,wzw1,wxu2,wyv2,wzw2
      integer :: s
      ! Two-stage RK2
      do s=1,2
         ! Pick interpolation point
         if (s.eq.1) then; pm=p1; else; pm=0.5_WP*(p1+p2); end if
         ! Cell-centered indices
         ipc=floor((pm(1)-xlo)*dxi-0.5_WP)
         jpc=floor((pm(2)-ylo)*dyi-0.5_WP)
         kpc=floor((pm(3)-zlo)*dzi-0.5_WP)
         ! Face-centered indices
         ipu=floor((pm(1)-xlo)*dxi)
         jpv=floor((pm(2)-ylo)*dyi)
         kpw=floor((pm(3)-zlo)*dzi)
         ! Clamp to array bounds
         ipu=max(lbound(pU,1),min(ubound(pU,1)-1,ipu))
         jpc=max(lbound(pU,2),min(ubound(pU,2)-1,jpc))
         kpc=max(lbound(pU,3),min(ubound(pU,3)-1,kpc))
         ipc=max(lbound(pV,1),min(ubound(pV,1)-1,ipc))
         jpv=max(lbound(pV,2),min(ubound(pV,2)-1,jpv))
         kpw=max(lbound(pW,3),min(ubound(pW,3)-1,kpw))
         ! Cell-centered weights
         wxc1=(pm(1)-(xlo+(real(ipc,WP)+0.5_WP)/dxi))*dxi; wxc1=max(0.0_WP,min(1.0_WP,wxc1)); wxc2=1.0_WP-wxc1
         wyc1=(pm(2)-(ylo+(real(jpc,WP)+0.5_WP)/dyi))*dyi; wyc1=max(0.0_WP,min(1.0_WP,wyc1)); wyc2=1.0_WP-wyc1
         wzc1=(pm(3)-(zlo+(real(kpc,WP)+0.5_WP)/dzi))*dzi; wzc1=max(0.0_WP,min(1.0_WP,wzc1)); wzc2=1.0_WP-wzc1
         ! Face-centered weights
         wxu1=(pm(1)-(xlo+real(ipu,WP)/dxi))*dxi; wxu1=max(0.0_WP,min(1.0_WP,wxu1)); wxu2=1.0_WP-wxu1
         wyv1=(pm(2)-(ylo+real(jpv,WP)/dyi))*dyi; wyv1=max(0.0_WP,min(1.0_WP,wyv1)); wyv2=1.0_WP-wyv1
         wzw1=(pm(3)-(zlo+real(kpw,WP)/dzi))*dzi; wzw1=max(0.0_WP,min(1.0_WP,wzw1)); wzw2=1.0_WP-wzw1
         ! U at x-faces: face-centered in x, cell-centered in y,z
         vel(1)=wzc1*(wyc1*(wxu1*pU(ipu+1,jpc+1,kpc+1,1)+wxu2*pU(ipu,jpc+1,kpc+1,1))+ &
         &            wyc2*(wxu1*pU(ipu+1,jpc  ,kpc+1,1)+wxu2*pU(ipu,jpc  ,kpc+1,1)))+&
         &      wzc2*(wyc1*(wxu1*pU(ipu+1,jpc+1,kpc  ,1)+wxu2*pU(ipu,jpc+1,kpc  ,1))+ &
         &            wyc2*(wxu1*pU(ipu+1,jpc  ,kpc  ,1)+wxu2*pU(ipu,jpc  ,kpc  ,1)))
         ! V at y-faces: cell-centered in x, face-centered in y, cell-centered in z
         vel(2)=wzc1*(wyv1*(wxc1*pV(ipc+1,jpv+1,kpc+1,1)+wxc2*pV(ipc,jpv+1,kpc+1,1))+ &
         &            wyv2*(wxc1*pV(ipc+1,jpv  ,kpc+1,1)+wxc2*pV(ipc,jpv  ,kpc+1,1)))+&
         &      wzc2*(wyv1*(wxc1*pV(ipc+1,jpv+1,kpc  ,1)+wxc2*pV(ipc,jpv+1,kpc  ,1))+ &
         &            wyv2*(wxc1*pV(ipc+1,jpv  ,kpc  ,1)+wxc2*pV(ipc,jpv  ,kpc  ,1)))
         ! W at z-faces: cell-centered in x,y, face-centered in z
         vel(3)=wzw1*(wyc1*(wxc1*pW(ipc+1,jpc+1,kpw+1,1)+wxc2*pW(ipc,jpc+1,kpw+1,1))+ &
         &            wyc2*(wxc1*pW(ipc+1,jpc  ,kpw+1,1)+wxc2*pW(ipc,jpc  ,kpw+1,1)))+&
         &      wzw2*(wyc1*(wxc1*pW(ipc+1,jpc+1,kpw  ,1)+wxc2*pW(ipc,jpc+1,kpw  ,1))+ &
         &            wyc2*(wxc1*pW(ipc+1,jpc  ,kpw  ,1)+wxc2*pW(ipc,jpc  ,kpw  ,1)))
         ! Project
         p2=p1+dt*vel
      end do
   end function project_rk2_staggered

   !> Cut a tet by PLIC plane and compute phasic flux (base case)
   !> Returns flux(8) = [Lvol, Gvol, Lbar(3), Gbar(3)]
   pure function tet2flux_plic(mytet,i0,j0,k0,pPLIC,plo) result(myflux)
      implicit none
      real(WP), dimension(3,4), intent(in) :: mytet
      integer, intent(in) :: i0,j0,k0
      integer, dimension(4), intent(in) :: plo
      real(WP), dimension(plo(1):,plo(2):,plo(3):,plo(4):), intent(in) :: pPLIC
      real(WP), dimension(8) :: myflux
      integer :: icase,n1,v1,v2
      real(WP), dimension(4) :: dd
      real(WP), dimension(3,8) :: vert
      real(WP), dimension(3) :: a,b,c,bary,normal
      real(WP) :: mu,my_vol,dist

      myflux=0.0_WP

      ! Pure cell shortcut - skip PLIC cutting
      if (pPLIC(i0,j0,k0,4).gt.+1.0e9_WP) then
         ! Pure liquid - all volume goes to liquid phase
         my_vol=abs(tet_vol(mytet))
         bary=0.25_WP*(mytet(:,1)+mytet(:,2)+mytet(:,3)+mytet(:,4))
         myflux(1)=my_vol
         myflux(3:5)=my_vol*bary
         return
      else if (pPLIC(i0,j0,k0,4).lt.-1.0e9_WP) then
         ! Pure gas - all volume goes to gas phase
         my_vol=abs(tet_vol(mytet))
         bary=0.25_WP*(mytet(:,1)+mytet(:,2)+mytet(:,3)+mytet(:,4))
         myflux(2)=my_vol
         myflux(6:8)=my_vol*bary
         return
      end if

      ! Get PLIC from this cell
      normal=pPLIC(i0,j0,k0,1:3)
      dist  =pPLIC(i0,j0,k0,4)

      ! Compute signed distance to plane for each vertex
      dd(1)=normal(1)*mytet(1,1)+normal(2)*mytet(2,1)+normal(3)*mytet(3,1)-dist
      dd(2)=normal(1)*mytet(1,2)+normal(2)*mytet(2,2)+normal(3)*mytet(3,2)-dist
      dd(3)=normal(1)*mytet(1,3)+normal(2)*mytet(2,3)+normal(3)*mytet(3,3)-dist
      dd(4)=normal(1)*mytet(1,4)+normal(2)*mytet(2,4)+normal(3)*mytet(3,4)-dist

      ! Find cut case
      icase=1+int(0.5_WP+sign(0.5_WP,dd(1)))+2*int(0.5_WP+sign(0.5_WP,dd(2))) &
      &    +4*int(0.5_WP+sign(0.5_WP,dd(3)))+8*int(0.5_WP+sign(0.5_WP,dd(4)))

      ! Copy vertices
      vert(:,1:4)=mytet(:,1:4)

      ! Create interpolated vertices on cut plane
      do n1=1,cut_nvert(icase)
         v1=cut_v1(n1,icase); v2=cut_v2(n1,icase)
         mu=min(1.0_WP,max(0.0_WP,-dd(v1)/(sign(abs(dd(v2)-dd(v1))+epsilon(1.0_WP),dd(v2)-dd(v1)))))
         vert(:,4+n1)=(1.0_WP-mu)*vert(:,v1)+mu*vert(:,v2)
      end do

      ! Gas tets: from 1 to cut_nntet-1
      do n1=1,cut_nntet(icase)-1
         a=vert(:,cut_vtet(1,n1,icase))-vert(:,cut_vtet(4,n1,icase))
         b=vert(:,cut_vtet(2,n1,icase))-vert(:,cut_vtet(4,n1,icase))
         c=vert(:,cut_vtet(3,n1,icase))-vert(:,cut_vtet(4,n1,icase))
         my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
         bary=0.25_WP*(vert(:,cut_vtet(1,n1,icase))+vert(:,cut_vtet(2,n1,icase)) &
         &            +vert(:,cut_vtet(3,n1,icase))+vert(:,cut_vtet(4,n1,icase)))
         myflux(2)=myflux(2)+my_vol
         myflux(6:8)=myflux(6:8)+my_vol*bary
      end do

      ! Liquid tets: from cut_ntets down to cut_nntet
      do n1=cut_ntets(icase),cut_nntet(icase),-1
         a=vert(:,cut_vtet(1,n1,icase))-vert(:,cut_vtet(4,n1,icase))
         b=vert(:,cut_vtet(2,n1,icase))-vert(:,cut_vtet(4,n1,icase))
         c=vert(:,cut_vtet(3,n1,icase))-vert(:,cut_vtet(4,n1,icase))
         my_vol=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
         bary=0.25_WP*(vert(:,cut_vtet(1,n1,icase))+vert(:,cut_vtet(2,n1,icase)) &
         &            +vert(:,cut_vtet(3,n1,icase))+vert(:,cut_vtet(4,n1,icase)))
         myflux(1)=myflux(1)+my_vol
         myflux(3:5)=myflux(3:5)+my_vol*bary
      end do

   end function tet2flux_plic

   !> Recursive grid-plane cutting: cut tet by grid planes, then by PLIC at leaf
   !> Returns flux(8) = [Lvol, Gvol, Lbar(3), Gbar(3)]
   pure recursive function tet2flux_recursive(mytet,myind,pPLIC,plo,xlo,ylo,zlo,dx,dy,dz) result(myflux)
      implicit none
      real(WP), dimension(3,4), intent(in) :: mytet
      integer,  dimension(3,4), intent(in) :: myind
      integer,  dimension(4),   intent(in) :: plo
      real(WP), dimension(plo(1):,plo(2):,plo(3):,plo(4):), intent(in) :: pPLIC
      real(WP), intent(in) :: xlo,ylo,zlo,dx,dy,dz
      real(WP), dimension(8) :: myflux
      real(WP) :: vol
      integer :: dir,cut_ind_,icase,n1,n2,v1,v2
      real(WP), dimension(4) :: dd
      real(WP), dimension(3,8) :: vert
      integer,  dimension(3,8,2) :: vert_ind
      real(WP) :: mu,my_vol_,xcut,ycut,zcut
      real(WP), dimension(3,4) :: newtet
      integer,  dimension(3,4) :: newind
      real(WP), dimension(3) :: a,b,c

      myflux=0.0_WP
      vol=dx*dy*dz

      ! Determine if tet spans multiple cells and needs cutting
      if (maxval(myind(1,:))-minval(myind(1,:)).gt.0) then
         dir=1; cut_ind_=maxval(myind(1,:))
         xcut=xlo+real(cut_ind_,WP)*dx
         dd(:)=mytet(1,:)-xcut
      else if (maxval(myind(2,:))-minval(myind(2,:)).gt.0) then
         dir=2; cut_ind_=maxval(myind(2,:))
         ycut=ylo+real(cut_ind_,WP)*dy
         dd(:)=mytet(2,:)-ycut
      else if (maxval(myind(3,:))-minval(myind(3,:)).gt.0) then
         dir=3; cut_ind_=maxval(myind(3,:))
         zcut=zlo+real(cut_ind_,WP)*dz
         dd(:)=mytet(3,:)-zcut
      else
         ! All vertices in same cell - cut by PLIC and return
         myflux=tet2flux_plic(mytet,myind(1,1),myind(2,1),myind(3,1),pPLIC,plo)
         return
      end if

      ! Find cut case (1-indexed: 1-16)
      icase=1+int(0.5_WP+sign(0.5_WP,dd(1)))+2*int(0.5_WP+sign(0.5_WP,dd(2))) &
      &    +4*int(0.5_WP+sign(0.5_WP,dd(3)))+8*int(0.5_WP+sign(0.5_WP,dd(4)))

      ! Copy vertices and indices
      do n1=1,4
         vert(:,n1)=mytet(:,n1)
         vert_ind(:,n1,1)=myind(:,n1)
         vert_ind(:,n1,2)=myind(:,n1)
         vert_ind(dir,n1,1)=min(vert_ind(dir,n1,1),cut_ind_-1)
         vert_ind(dir,n1,2)=max(vert_ind(dir,n1,1),cut_ind_)
      end do

      ! Create interpolated vertices on cut plane
      do n1=1,cut_nvert(icase)
         v1=cut_v1(n1,icase); v2=cut_v2(n1,icase)
         mu=min(1.0_WP,max(0.0_WP,-dd(v1)/(sign(abs(dd(v2)-dd(v1))+epsilon(1.0_WP),dd(v2)-dd(v1)))))
         vert(:,4+n1)=(1.0_WP-mu)*vert(:,v1)+mu*vert(:,v2)
         vert_ind(1,4+n1,1)=floor((vert(1,4+n1)-xlo)/dx)
         vert_ind(2,4+n1,1)=floor((vert(2,4+n1)-ylo)/dy)
         vert_ind(3,4+n1,1)=floor((vert(3,4+n1)-zlo)/dz)
         vert_ind(:,4+n1,1)=max(vert_ind(:,4+n1,1),min(vert_ind(:,v1,1),vert_ind(:,v2,1)))
         vert_ind(:,4+n1,1)=min(vert_ind(:,4+n1,1),max(vert_ind(:,v1,1),vert_ind(:,v2,1)))
         vert_ind(:,4+n1,2)=vert_ind(:,4+n1,1)
         vert_ind(dir,4+n1,1)=cut_ind_-1
         vert_ind(dir,4+n1,2)=cut_ind_
      end do

      ! Create and process sub-tets
      do n1=1,cut_ntets(icase)
         do n2=1,4
            newtet(:,n2)=vert(:,cut_vtet(n2,n1,icase))
            newind(:,n2)=vert_ind(:,cut_vtet(n2,n1,icase),cut_side(n1,icase))
         end do
         ! Check for zero-volume tet
         a=newtet(:,1)-newtet(:,4)
         b=newtet(:,2)-newtet(:,4)
         c=newtet(:,3)-newtet(:,4)
         my_vol_=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
         if (my_vol_.lt.1.0e-15_WP*vol) cycle
         myflux=myflux+tet2flux_recursive(newtet,newind,pPLIC,plo,xlo,ylo,zlo,dx,dy,dz)
      end do

   end function tet2flux_recursive

   !> Iterative (stack-based) grid-plane cutting: same result as recursive version
   !> GPU-friendly: no recursion, uses explicit stack
   !> Returns flux(8) = [Lvol, Gvol, Lbar(3), Gbar(3)]
   pure subroutine tet2flux_iterative(mytet,myind,pPLIC,plo,xlo,ylo,zlo,dx,dy,dz,myflux)
      implicit none
      real(WP), dimension(3,4), intent(in) :: mytet
      integer,  dimension(3,4), intent(in) :: myind
      integer,  dimension(4),   intent(in) :: plo
      real(WP), dimension(plo(1):,plo(2):,plo(3):,plo(4):), intent(in) :: pPLIC
      real(WP), intent(in) :: xlo,ylo,zlo,dx,dy,dz
      real(WP), dimension(8), intent(out) :: myflux
      ! Stack parameters
      integer, parameter :: STACK_MAX=64
      real(WP), dimension(3,4,STACK_MAX) :: stack_tet
      integer,  dimension(3,4,STACK_MAX) :: stack_ind
      integer :: sp
      ! Working variables
      real(WP), dimension(3,4) :: cur_tet,newtet
      integer,  dimension(3,4) :: cur_ind,newind
      integer  :: dir,cut_ind_,icase,n1,n2,v1,v2
      real(WP), dimension(4) :: dd
      real(WP), dimension(3,8) :: vert
      integer,  dimension(3,8,2) :: vert_ind
      real(WP) :: mu,my_vol_,vol
      real(WP), dimension(3) :: a,b,c

      myflux=0.0_WP
      vol=dx*dy*dz

      ! Push initial tet onto stack
      sp=1
      stack_tet(:,:,1)=mytet
      stack_ind(:,:,1)=myind

      ! Process stack
      do while (sp.gt.0)
         ! Pop
         cur_tet=stack_tet(:,:,sp)
         cur_ind=stack_ind(:,:,sp)
         sp=sp-1

         ! Determine if tet spans multiple cells
         if (maxval(cur_ind(1,:))-minval(cur_ind(1,:)).gt.0) then
            dir=1; cut_ind_=maxval(cur_ind(1,:))
            dd(:)=cur_tet(1,:)-(xlo+real(cut_ind_,WP)*dx)
         else if (maxval(cur_ind(2,:))-minval(cur_ind(2,:)).gt.0) then
            dir=2; cut_ind_=maxval(cur_ind(2,:))
            dd(:)=cur_tet(2,:)-(ylo+real(cut_ind_,WP)*dy)
         else if (maxval(cur_ind(3,:))-minval(cur_ind(3,:)).gt.0) then
            dir=3; cut_ind_=maxval(cur_ind(3,:))
            dd(:)=cur_tet(3,:)-(zlo+real(cut_ind_,WP)*dz)
         else
            ! Leaf: cut by PLIC
            myflux=myflux+tet2flux_plic(cur_tet,cur_ind(1,1),cur_ind(2,1),cur_ind(3,1),pPLIC,plo)
            cycle
         end if

         ! Find cut case
         icase=1+int(0.5_WP+sign(0.5_WP,dd(1)))+2*int(0.5_WP+sign(0.5_WP,dd(2))) &
         &    +4*int(0.5_WP+sign(0.5_WP,dd(3)))+8*int(0.5_WP+sign(0.5_WP,dd(4)))

         ! Copy vertices and indices
         do n1=1,4
            vert(:,n1)=cur_tet(:,n1)
            vert_ind(:,n1,1)=cur_ind(:,n1)
            vert_ind(:,n1,2)=cur_ind(:,n1)
            vert_ind(dir,n1,1)=min(vert_ind(dir,n1,1),cut_ind_-1)
            vert_ind(dir,n1,2)=max(vert_ind(dir,n1,1),cut_ind_)
         end do

         ! Create interpolated vertices on cut plane
         do n1=1,cut_nvert(icase)
            v1=cut_v1(n1,icase); v2=cut_v2(n1,icase)
            mu=min(1.0_WP,max(0.0_WP,-dd(v1)/(sign(abs(dd(v2)-dd(v1))+epsilon(1.0_WP),dd(v2)-dd(v1)))))
            vert(:,4+n1)=(1.0_WP-mu)*vert(:,v1)+mu*vert(:,v2)
            vert_ind(1,4+n1,1)=floor((vert(1,4+n1)-xlo)/dx)
            vert_ind(2,4+n1,1)=floor((vert(2,4+n1)-ylo)/dy)
            vert_ind(3,4+n1,1)=floor((vert(3,4+n1)-zlo)/dz)
            vert_ind(:,4+n1,1)=max(vert_ind(:,4+n1,1),min(vert_ind(:,v1,1),vert_ind(:,v2,1)))
            vert_ind(:,4+n1,1)=min(vert_ind(:,4+n1,1),max(vert_ind(:,v1,1),vert_ind(:,v2,1)))
            vert_ind(:,4+n1,2)=vert_ind(:,4+n1,1)
            vert_ind(dir,4+n1,1)=cut_ind_-1
            vert_ind(dir,4+n1,2)=cut_ind_
         end do

         ! Create sub-tets and push onto stack
         do n1=1,cut_ntets(icase)
            do n2=1,4
               newtet(:,n2)=vert(:,cut_vtet(n2,n1,icase))
               newind(:,n2)=vert_ind(:,cut_vtet(n2,n1,icase),cut_side(n1,icase))
            end do
            ! Check for zero-volume tet
            a=newtet(:,1)-newtet(:,4)
            b=newtet(:,2)-newtet(:,4)
            c=newtet(:,3)-newtet(:,4)
            my_vol_=abs(a(1)*(b(2)*c(3)-c(2)*b(3))-a(2)*(b(1)*c(3)-c(1)*b(3))+a(3)*(b(1)*c(2)-c(1)*b(2)))/6.0_WP
            if (my_vol_.lt.1.0e-15_WP*vol) cycle
            ! Push onto stack
            sp=sp+1
            if (sp.gt.STACK_MAX) return ! Safety: bail if stack overflows
            stack_tet(:,:,sp)=newtet
            stack_ind(:,:,sp)=newind
         end do
      end do

   end subroutine tet2flux_iterative

end module amrvof_geometry