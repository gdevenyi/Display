
#include  <display.h>

private  void  intersect_plane_polygons(
    Vector            *plane_normal,
    Real              plane_constant,
    polygons_struct   *polygons,
    lines_struct      *lines,
    int               *n_points_alloced,
    int               *n_indices_alloced,
    int               *n_end_indices_alloced );

public  void  intersect_plane_with_polygons(
    display_struct    *display,
    Vector            *plane_normal,
    Real              plane_constant,
    lines_struct      *lines,
    int               *n_points_alloced,
    int               *n_indices_alloced,
    int               *n_end_indices_alloced )
{
    object_struct            *object;
    object_traverse_struct   object_traverse;

    lines->n_items = 0;
    lines->n_points = 0;

    initialize_object_traverse( &object_traverse, N_MODELS, display->models );

    while( get_next_object_traverse(&object_traverse,&object) )
    {
        if( object->object_type == POLYGONS && object->visibility )
        {
            intersect_plane_polygons( plane_normal, plane_constant,
                                      get_polygons_ptr(object), lines,
                                      n_points_alloced, n_indices_alloced,
                                      n_end_indices_alloced );
        }
    }
}

private  void  intersect_plane_polygons(
    Vector            *plane_normal,
    Real              plane_constant,
    polygons_struct   *polygons,
    lines_struct      *lines,
    int               *n_points_alloced,
    int               *n_indices_alloced,
    int               *n_end_indices_alloced )
{
    int       i;

    if( FALSE && polygons->bintree != (bintree_struct *) 0 )
    {
/*
        intersect_ray_with_bintree( ray_origin, ray_direction,
                                                 polygons->bintree, polygons,
                                                 poly_index, dist );
*/
    }
    else
    {
        for_less( i, 0, polygons->n_items )
        {
            (void) intersect_plane_one_polygon( plane_normal, plane_constant,
                                                polygons, i, lines, 
                                                n_points_alloced,
                                                n_indices_alloced,
                                                n_end_indices_alloced );
        }
    }
}

public  BOOLEAN  intersect_plane_one_polygon(
    Vector            *plane_normal,
    Real              plane_constant,
    polygons_struct   *polygons,
    int               poly,
    lines_struct      *lines,
    int               *n_points_alloced,
    int               *n_indices_alloced,
    int               *n_end_indices_alloced )
{
    int       point_index, n_indices;
    Point     points[2];
    BOOLEAN   intersects;

    intersects = get_plane_polygon_intersection( plane_normal, plane_constant,
                                                 polygons, poly, points );

    if( intersects )
    {
        point_index = lines->n_points;

        ADD_ELEMENT_TO_ARRAY_WITH_SIZE( lines->points,
                                        *n_points_alloced, lines->n_points,
                                        points[0], DEFAULT_CHUNK_SIZE );
        ADD_ELEMENT_TO_ARRAY_WITH_SIZE( lines->points,
                                        *n_points_alloced, lines->n_points,
                                        points[1], DEFAULT_CHUNK_SIZE );

        n_indices = NUMBER_INDICES( *lines );

        ADD_ELEMENT_TO_ARRAY_WITH_SIZE( lines->indices,
                                        *n_indices_alloced, n_indices,
                                        point_index, DEFAULT_CHUNK_SIZE );
        ADD_ELEMENT_TO_ARRAY_WITH_SIZE( lines->indices,
                                        *n_indices_alloced, n_indices,
                                        point_index+1, DEFAULT_CHUNK_SIZE );

        ADD_ELEMENT_TO_ARRAY_WITH_SIZE( lines->end_indices,
                                        *n_end_indices_alloced, lines->n_items,
                                        n_indices, DEFAULT_CHUNK_SIZE );
    }

    return( intersects );
}
