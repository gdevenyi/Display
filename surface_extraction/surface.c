/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/visualization/Display/surface_extraction/surface.c,v 1.53 1995-10-19 15:52:38 david Exp $";
#endif


#include  <display.h>

private  BOOLEAN  find_close_voxel_containing_range(
    Volume                     volume,
    Volume                     label_volume,
    unsigned_byte              voxel_done_flags[],
    surface_extraction_struct  *surface_extraction,
    int                        x,
    int                        y,
    int                        z,
    voxel_index_struct         *found_indices );
private  void  add_voxel_neighbours(
    Volume                              volume,
    Volume                              label_volume,
    int                                 sizes[],
    int                                 x,
    int                                 y,
    int                                 z,
    BOOLEAN                             surface_only,
    surface_extraction_struct           *surface_extraction,
    bitlist_3d_struct                   *voxels_queued,
    voxel_queue_struct                  *voxel_queue );
private  BOOLEAN  cube_is_within_distance(
    surface_extraction_struct           *surface_extraction,
    int                                 x,
    int                                 y,
    int                                 z );
private  void  possibly_output(
    polygons_struct  *p );
private  void  delete_edge_points_no_longer_needed(
    Volume              volume,
    voxel_index_struct  *voxel_index,
    unsigned_byte       voxel_done_flags[],
    hash_table_struct   *edge_points );

private  BOOLEAN  surface_voxel_is_within_volume(
    int       sizes[],
    int       indices[] )
{
    return( indices[X] >= 0 && indices[X] < sizes[X]-1 &&
            indices[Y] >= 0 && indices[Y] < sizes[Y]-1 &&
            indices[Z] >= 0 && indices[Z] < sizes[Z]-1 );
}

public  void  start_surface_extraction_at_point(
    display_struct     *display,
    Volume             volume,
    Volume             label_volume,
    BOOLEAN            binary_flag,
    Real               min_value,
    Real               max_value,
    int                x,
    int                y,
    int                z )
{
    int                         sizes[N_DIMENSIONS];
    int                         indices[N_DIMENSIONS];
    surface_extraction_struct   *surface_extraction;
    voxel_index_struct          voxel_indices;

    initialize_surface_extraction_for_volume( display, volume, label_volume );

    surface_extraction = &display->three_d.surface_extraction;
    surface_extraction->binary_flag = binary_flag;
    surface_extraction->min_value = min_value;
    surface_extraction->max_value = max_value;

    indices[X] = x;
    indices[Y] = y;
    indices[Z] = z;
    get_volume_sizes( volume, sizes );
    if( surface_voxel_is_within_volume( sizes, indices ) )
    {
        display->three_d.surface_extraction.x_starting_voxel = x;
        display->three_d.surface_extraction.y_starting_voxel = y;
        display->three_d.surface_extraction.z_starting_voxel = z;

        while( voxels_remaining(
                 &display->three_d.surface_extraction.voxels_to_do ) )
        {
            get_next_voxel_from_queue( 
                 &display->three_d.surface_extraction.voxels_to_do,
                 &voxel_indices );

            reset_voxel_flag( 
                     &display->three_d.surface_extraction.voxels_queued,
                     &voxel_indices );
        }

        if( find_close_voxel_containing_range( volume, label_volume,
                  display->three_d.surface_extraction.voxel_done_flags,
                  &display->three_d.surface_extraction,
                  x, y, z, &voxel_indices ) )
        {
            insert_in_voxel_queue(
                         &display->three_d.surface_extraction.voxels_to_do,
                         &voxel_indices );

            set_voxel_flag( &display->three_d.surface_extraction.voxels_queued,
                            &voxel_indices );

            start_surface_extraction( display );
        }
        else
            print( "No surface voxels found.\n" );
    }
}

private  BOOLEAN  find_close_voxel_containing_range(
    Volume                     volume,
    Volume                     label_volume,
    unsigned_byte              voxel_done_flags[],
    surface_extraction_struct  *surface_extraction,
    int                        x,
    int                        y,
    int                        z,
    voxel_index_struct         *found_indices )
{
    BOOLEAN                   found, voxel_contains;
    int                       sizes[MAX_DIMENSIONS], voxel[MAX_DIMENSIONS];
    BOOLEAN                   voxel_done;
    voxel_queue_struct        voxels_to_check;
    voxel_index_struct        indices, insert;
    bitlist_3d_struct         voxels_searched;

    get_volume_sizes( volume, sizes );

    insert.i[X] = (short) MIN( x, sizes[X]-2 );
    insert.i[Y] = (short) MIN( y, sizes[Y]-2 );
    insert.i[Z] = (short) MIN( z, sizes[Z]-2 );

    found = FALSE;

    initialize_voxel_flags( &voxels_searched, sizes );

    initialize_voxel_queue( &voxels_to_check );

    insert_in_voxel_queue( &voxels_to_check, &insert );

    set_voxel_flag( &voxels_searched, &insert );

    while( !found && voxels_remaining(&voxels_to_check) )
    {
        get_next_voxel_from_queue( &voxels_to_check, &indices );

        voxel[X] = indices.i[X];
        voxel[Y] = indices.i[Y];
        voxel[Z] = indices.i[Z];
        voxel_contains = voxel_contains_surface( volume, label_volume,
                                                 surface_extraction, voxel );

        voxel_done = get_voxel_done_flag( sizes, voxel_done_flags, &indices );

        if( voxel_contains && voxel_done == 0 )
        {
            found_indices->i[X] = indices.i[X];
            found_indices->i[Y] = indices.i[Y];
            found_indices->i[Z] = indices.i[Z];
            found = TRUE;
        }
        else if( voxel_contains || voxel_done == 0 )
        {
            add_voxel_neighbours( volume, label_volume, sizes,
                                  indices.i[X], indices.i[Y], indices.i[Z],
                                  voxel_done, surface_extraction,
                                  &voxels_searched, &voxels_to_check );
        }
    }

    delete_voxel_flags( &voxels_searched );

    delete_voxel_queue( &voxels_to_check );

    return( found );
}

public  void  extract_more_surface(
    display_struct    *display )
{
    int                         n_voxels_done, sizes[N_DIMENSIONS];
    voxel_index_struct          voxel_index;
    surface_extraction_struct   *surface_extraction;
    Volume                      volume, label_volume;
    Real                        stop_time;

    n_voxels_done = 0;

    surface_extraction = &display->three_d.surface_extraction;
    volume = surface_extraction->volume;
    label_volume = surface_extraction->label_volume;

    get_volume_sizes( volume, sizes );

    stop_time = current_realtime_seconds() + Max_seconds_per_voxel_update;

    while( (n_voxels_done < Min_voxels_per_update ||
           (n_voxels_done < Max_voxels_per_update &&
            current_realtime_seconds() < stop_time) ) &&
           voxels_remaining( &surface_extraction->voxels_to_do ) )
    {
        possibly_output( surface_extraction->polygons );

        get_next_voxel_from_queue( &surface_extraction->voxels_to_do,
                                   &voxel_index );

        reset_voxel_flag( &surface_extraction->voxels_queued, &voxel_index);

        if( extract_voxel_surface( volume, label_volume,
                                   surface_extraction, &voxel_index,
                            surface_extraction->n_voxels_with_surface == 0) )
        {
            ++n_voxels_done;
            ++surface_extraction->n_voxels_with_surface;

            delete_edge_points_no_longer_needed( volume, &voxel_index,
                                 surface_extraction->voxel_done_flags,
                                 &surface_extraction->edge_points );

            add_voxel_neighbours( volume, label_volume, sizes,
                        voxel_index.i[X], voxel_index.i[Y], voxel_index.i[Z],
                        TRUE,
                        surface_extraction,
                        &surface_extraction->voxels_queued,
                        &surface_extraction->voxels_to_do );
        }
    }

    if( !voxels_remaining( &surface_extraction->voxels_to_do ) )
    {
        print( "Surface extraction finished\n" );
        stop_surface_extraction( display );
        update_all_menu_text( display );
    }
}

private  void  add_voxel_neighbours(
    Volume                          volume,
    Volume                          label_volume,
    int                             sizes[],
    int                             x,
    int                             y,
    int                             z,
    BOOLEAN                         surface_only,
    surface_extraction_struct       *surface_extraction,
    bitlist_3d_struct               *voxels_queued,
    voxel_queue_struct              *voxel_queue )
{
    int                      x_offset, y_offset, z_offset;
    int                      indices[N_DIMENSIONS];
    voxel_index_struct       neighbour;

    for_inclusive( x_offset, -1, 1 )
    {
        neighbour.i[X] = x + x_offset;

        for_inclusive( y_offset, -1, 1 )
        {
            neighbour.i[Y] = y + y_offset;
            for_inclusive( z_offset, -1, 1 )
            {
                neighbour.i[Z] = z + z_offset;
                indices[X] = neighbour.i[X];
                indices[Y] = neighbour.i[Y];
                indices[Z] = neighbour.i[Z];
                if( (x_offset != 0 || y_offset != 0 || z_offset != 0) &&
                    surface_voxel_is_within_volume( sizes, indices ) &&
                    !get_voxel_flag( voxels_queued, &neighbour ) &&
                    get_voxel_done_flag( sizes,
                                         surface_extraction->voxel_done_flags,
                                        &neighbour ) != VOXEL_COMPLETELY_DONE &&
                    cube_is_within_distance( surface_extraction,
                                     indices[X], indices[Y], indices[Z] ) )
                {
                    set_voxel_flag( voxels_queued, &neighbour );
                    if( !surface_only ||
                        voxel_contains_surface( volume, label_volume,
                                                surface_extraction, indices ) )
                    {
                        insert_in_voxel_queue( voxel_queue, &neighbour );
                    }
                }
            }
        }
    }
}

private  BOOLEAN  cube_is_within_distance(
    surface_extraction_struct           *surface_extraction,
    int                                 x,
    int                                 y,
    int                                 z )
{
    BOOLEAN  within_dist;
    int      dx, dy, dz;

    dx = ABS( x - surface_extraction->x_starting_voxel );
    dy = ABS( y - surface_extraction->y_starting_voxel );
    dz = ABS( z - surface_extraction->z_starting_voxel );

    within_dist = dx <= surface_extraction->x_voxel_max_distance &&
                  dy <= surface_extraction->y_voxel_max_distance &&
                  dz <= surface_extraction->z_voxel_max_distance;

    return( within_dist );
}

private  void  possibly_output(
    polygons_struct  *p )
{
    static   int  count  = 0;
    Status   status;
    FILE     *file;
    char     buffer[EXTREMELY_LARGE_STRING_SIZE];

    if( Output_every > 0 )
    {
        if( (count-1) * Output_every > p->n_items )
        {
            count = p->n_items / Output_every;
        }

        if( p->n_items > count * Output_every )
        {
            ++count;
            (void) sprintf( buffer, Tmp_surface_name, count );

            status = open_file_with_default_suffix( buffer, "obj",
                                            WRITE_FILE, BINARY_FORMAT, &file );

            if( status == OK )
            {
                status = io_polygons( file, WRITE_FILE, BINARY_FORMAT, p );
            }

            if( status == OK )
            {
                status = close_file( file );
            }
        }
    }
}

private  void  delete_edge_points_no_longer_needed(
    Volume              volume,
    voxel_index_struct  *voxel_index,
    unsigned_byte       voxel_done_flags[],
    hash_table_struct   *edge_points )
{
    int                 axis_index, a1, a2;
    int                 x, y, dx, dy, dz;
    BOOLEAN             all_four_done;
    BOOLEAN             voxel_done[3][3][3];
    voxel_index_struct  indices;
    int                 int_indices[N_DIMENSIONS];
    int                 sizes[N_DIMENSIONS];

    get_volume_sizes( volume, sizes );

    for_inclusive( dx, -1, 1 )
    {
        indices.i[X] = voxel_index->i[X] + dx;
        for_inclusive( dy, -1, 1 )
        {
            indices.i[Y] = voxel_index->i[Y] + dy;
            for_inclusive( dz, -1, 1 )
            {
                indices.i[Z] = voxel_index->i[Z] + dz;
                int_indices[X] = indices.i[X];
                int_indices[Y] = indices.i[Y];
                int_indices[Z] = indices.i[Z];

                if( !surface_voxel_is_within_volume( sizes, int_indices ) ||
                    get_voxel_done_flag( sizes, voxel_done_flags, &indices )
                    == VOXEL_COMPLETELY_DONE )
                {
                    voxel_done[dx+1][dy+1][dz+1] = TRUE;
                }
                else
                {
                    voxel_done[dx+1][dy+1][dz+1] = FALSE;
                }
            }
        }
    }

    for_less( axis_index, 0, N_DIMENSIONS )
    {
        a1 = (axis_index + 1) % N_DIMENSIONS;
        a2 = (axis_index + 2) % N_DIMENSIONS;

        for_less( x, 0, 2 )
        {
            for_less( y, 0, 2 )
            {
                all_four_done = TRUE;

                for_less( dx, 0, 2 )
                {
                    for_less( dy, 0, 2 )
                    {
                        indices.i[axis_index] = 1;
                        indices.i[a1] = x+dx;
                        indices.i[a2] = y+dy;

                        if( !voxel_done[indices.i[X]]
                                       [indices.i[Y]]
                                       [indices.i[Z]] )
                        {
                            all_four_done = FALSE;
                            break;
                        }
                    }
                }

                if( all_four_done )
                {
                    indices.i[axis_index] = voxel_index->i[axis_index];
                    indices.i[a1] = voxel_index->i[a1] + x;
                    indices.i[a2] = voxel_index->i[a2] + y;

                    remove_edge_point( sizes, edge_points, &indices,
                                       axis_index );
                }
            }
        }
    }
}
