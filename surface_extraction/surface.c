
#include  <def_graphics.h>
#include  <def_globals.h>
#include  <def_marching_cubes.h>
#include  <def_splines.h>
#include  <def_bitlist.h>
#include  <def_files.h>

static    Status   delete_edge_points_no_longer_needed();
static    void     possibly_output();

public  void  start_surface_extraction_at_point( graphics, x, y, z )
    graphics_struct    *graphics;
    int                x, y, z;
{
    surface_extraction_struct   *surface_extraction;
    Boolean                     cube_is_within_volume();
    Boolean                     find_close_voxel_containing_value();
    voxel_index_struct          voxel_indices;
    Status                      status;
    Status                      insert_in_voxel_queue();
    Status                      set_voxel_flag();
    Status                      reset_voxel_flag();
    void                        start_surface_extraction();
    void                        get_next_voxel_from_queue();
    void                        set_isosurface_value();

    surface_extraction = &graphics->three_d.surface_extraction;

    status = OK;

    if( cube_is_within_volume(
             graphics->associated[SLICE_WINDOW]->slice.volume, x, y, z ) )
    {
        graphics->three_d.surface_extraction.x_starting_voxel = x;
        graphics->three_d.surface_extraction.y_starting_voxel = y;
        graphics->three_d.surface_extraction.z_starting_voxel = z;

        if( !surface_extraction->isovalue_selected )
        {
            set_isosurface_value( surface_extraction );
        }
        else
        {
            while( voxels_remaining(
                     &graphics->three_d.surface_extraction.voxels_to_do ) )
            {
                get_next_voxel_from_queue( 
                     &graphics->three_d.surface_extraction.voxels_to_do,
                     &voxel_indices );

                if( status == OK )
                {
                    status = reset_voxel_flag(
                         graphics->associated[SLICE_WINDOW]->slice.volume,
                         &graphics->three_d.surface_extraction.voxels_queued,
                         &voxel_indices );
                }
            }
        }

        if( find_close_voxel_containing_value(
                  graphics->associated[SLICE_WINDOW]->slice.volume,
                  graphics->three_d.surface_extraction.voxel_done_flags,
                  graphics->three_d.surface_extraction.isovalue,
                  &graphics->three_d.surface_extraction,
                  x, y, z, &voxel_indices ) )
        {
            if( status == OK )
            {
                status = insert_in_voxel_queue(
                         &graphics->three_d.surface_extraction.voxels_to_do,
                         &voxel_indices );
            }

            if( status == OK )
            {
                status = set_voxel_flag( 
                           graphics->associated[SLICE_WINDOW]->slice.volume,
                           &graphics->three_d.surface_extraction.voxels_queued,
                           &voxel_indices );
            }

            start_surface_extraction( graphics );
        }
    }
}

private  Boolean  find_close_voxel_containing_value( volume, voxel_done_flags,
                                                    value, surface_extraction,
                                                    x, y, z, found_indices )
    volume_struct              *volume;
    unsigned_byte              voxel_done_flags[];
    Real                       value;
    surface_extraction_struct  *surface_extraction;
    int                        x, y, z;
    voxel_index_struct         *found_indices;
{
    Status                                status;
    Boolean                               found, voxel_contains;
    int                                   nx, ny, nz;
    unsigned_byte                         voxel_done;
    unsigned_byte                         get_voxel_done_flag();
    Boolean                               voxel_contains_value();
    QUEUE_STRUCT( voxel_index_struct )    voxels_to_check;
    voxel_index_struct                    indices, insert;
    bitlist_struct                        voxels_searched;
    Status                                set_voxel_flag();
    Status                                delete_voxel_flags();
    Status                                insert_in_voxel_queue();
    Status                                delete_voxel_queue();
    Status                                initialize_voxel_flags();
    void                                  initialize_voxel_queue();
    void                                  get_next_voxel_from_queue();
    void                                  add_voxel_neighbours();
    void                                  get_volume_size();

    get_volume_size( volume, &nx, &ny, &nz );

    insert.i[X] = MIN( x, nx-2 );
    insert.i[Y] = MIN( y, ny-2 );
    insert.i[Z] = MIN( z, nz-2 );

    found = FALSE;

    status = initialize_voxel_flags( &voxels_searched, get_n_voxels(volume) );

    initialize_voxel_queue( &voxels_to_check );

    if( status == OK )
    {
        status = insert_in_voxel_queue( &voxels_to_check, &insert );
    }

    if( status == OK )
    {
        status = set_voxel_flag( volume, &voxels_searched, &insert );
    }

    while( !found && status == OK && voxels_remaining(&voxels_to_check) )
    {
        get_next_voxel_from_queue( &voxels_to_check, &indices );

        voxel_contains = voxel_contains_value( volume,
                                               indices.i[X],
                                               indices.i[Y],
                                               indices.i[Z], value );

        voxel_done = get_voxel_done_flag( volume, voxel_done_flags, &indices );

        if( voxel_contains && !voxel_done )
        {
            found_indices->i[X] = indices.i[X];
            found_indices->i[Y] = indices.i[Y];
            found_indices->i[Z] = indices.i[Z];
            found = TRUE;
        }
        else if( voxel_contains || !voxel_done )
        {
            add_voxel_neighbours( volume,
                                  indices.i[X],
                                  indices.i[Y],
                                  indices.i[Z],
                                  (Boolean) voxel_done, value,
                                  surface_extraction,
                                  &voxels_searched, &voxels_to_check );
        }
    }

    if( status == OK )
    {
        status = delete_voxel_flags( &voxels_searched );
    }

    if( status == OK )
    {
        status = delete_voxel_queue( &voxels_to_check );
    }

    return( found );
}

public  Status  extract_more_surface( graphics )
    graphics_struct   *graphics;
{
    int                         n_voxels_done;
    voxel_index_struct          voxel_index;
    surface_extraction_struct   *surface_extraction;
    volume_struct               *volume;
    void                        get_next_voxel_from_queue();
    Boolean                     extract_voxel_surface();
    void                        add_voxel_neighbours();
    Real                        stop_time;
    Real                        current_realtime_seconds();
    Status                      status;
    Status                      set_voxel_flag();
    Status                      reset_voxel_flag();
    Status                      set_voxel_done_flag();
    void                        stop_surface_extraction();
    void                        label_voxel_as_done();

    status = OK;

    n_voxels_done = 0;

    surface_extraction = &graphics->three_d.surface_extraction;
    volume = graphics->associated[SLICE_WINDOW]->slice.volume;

    stop_time = current_realtime_seconds() + Max_seconds_per_voxel_update;

    while( status == OK && (n_voxels_done < Min_voxels_per_update ||
           (n_voxels_done < Max_voxels_per_update &&
            current_realtime_seconds() < stop_time) ) &&
           voxels_remaining( &surface_extraction->voxels_to_do ) )
    {
        possibly_output( surface_extraction->polygons );

        get_next_voxel_from_queue( &surface_extraction->voxels_to_do,
                                   &voxel_index );

        if( status == OK )
        {
            status = reset_voxel_flag( volume,
                                       &surface_extraction->voxels_queued,
                                       &voxel_index);
        }

        if( extract_voxel_surface( volume, surface_extraction, &voxel_index,
                            surface_extraction->n_voxels_with_surface == 0) )
        {
            ++n_voxels_done;
            ++surface_extraction->n_voxels_with_surface;

            if( status == OK )
            {
                status = delete_edge_points_no_longer_needed( 
                             volume, &voxel_index,
                             surface_extraction->voxel_done_flags,
                             &surface_extraction->edge_points );
            }

            if( status == OK && Display_surface_in_slices )
            {
                label_voxel_as_done( volume,
                                     voxel_index.i[X],
                                     voxel_index.i[Y],
                                     voxel_index.i[Z] );
            }

            if( status == OK )
            {
                add_voxel_neighbours(
                        volume,
                        voxel_index.i[X],
                        voxel_index.i[Y],
                        voxel_index.i[Z],
                        TRUE, surface_extraction->isovalue,
                        surface_extraction,
                        &surface_extraction->voxels_queued,
                        &surface_extraction->voxels_to_do );
            }
        }
    }

    if( !voxels_remaining( &surface_extraction->voxels_to_do ) )
    {
        PRINT( "Surface extraction finished\n" );
        stop_surface_extraction( graphics );
    }

    return( status );
}

private  void  add_voxel_neighbours( volume, x, y, z, surface_only, isovalue,
                                     surface_extraction,
                                     voxels_queued, voxel_queue )
    volume_struct                       *volume;
    int                                 x, y, z;
    Boolean                             surface_only;
    Real                                isovalue;
    surface_extraction_struct           *surface_extraction;
    bitlist_struct                      *voxels_queued;
    QUEUE_STRUCT(voxel_index_struct)    *voxel_queue;
{
    Status                   status;
    int                      x_offset, y_offset, z_offset;
    voxel_index_struct       neighbour;
    Boolean                  cube_is_within_volume();
    Status                   insert_in_voxel_queue();
    Status                   set_voxel_flag();

    status = OK;

    for_inclusive( x_offset, -1, 1 )
    {
        neighbour.i[X] = x + x_offset;

        for_inclusive( y_offset, -1, 1 )
        {
            neighbour.i[Y] = y + y_offset;
            for_inclusive( z_offset, -1, 1 )
            {
                neighbour.i[Z] = z + z_offset;
                if( (x_offset != 0 || y_offset != 0 || z_offset != 0) &&
                    cube_is_within_volume( volume,
                                           neighbour.i[X],
                                           neighbour.i[Y],
                                           neighbour.i[Z] ) &&
                    cube_is_within_distance( surface_extraction,
                                             neighbour.i[X],
                                             neighbour.i[Y],
                                             neighbour.i[Z] ) &&
                    !get_voxel_flag( volume, voxels_queued, &neighbour ) )
                {
                    status = set_voxel_flag( volume, voxels_queued, &neighbour);
                    if( status == OK &&
                        (!surface_only ||
                         voxel_contains_value( volume,
                                               neighbour.i[X],
                                               neighbour.i[Y],
                                               neighbour.i[Z],
                                               isovalue )) )
                    {
                        status = insert_in_voxel_queue( voxel_queue,&neighbour);
                    }
                }
            }
        }
    }
}

private  Boolean  cube_is_within_distance( surface_extraction, x, y, z )
    surface_extraction_struct           *surface_extraction;
    int                                 x, y, z;
{
    Boolean  within_dist;
    int      dx, dy, dz;

    dx = ABS( x - surface_extraction->x_starting_voxel );
    dy = ABS( y - surface_extraction->y_starting_voxel );
    dz = ABS( z - surface_extraction->z_starting_voxel );

    within_dist = dx <= surface_extraction->x_voxel_max_distance &&
                  dy <= surface_extraction->y_voxel_max_distance &&
                  dz <= surface_extraction->z_voxel_max_distance;

    return( within_dist );
}

private  void  possibly_output( p )
    polygons_struct  *p;
{
    static   int  count  = 0;
    Status   status;
    FILE     *file;
    String   name;
    Status   io_polygons();

    if( Output_every > 0 )
    {
        if( (count-1) * Output_every > p->n_items )
        {
            count = p->n_items / Output_every;
        }

        if( p->n_items > count * Output_every )
        {
            ++count;
            (void) sprintf( name, Tmp_surface_name, count );

            status = open_file_with_default_suffix( name, "obj",
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

private  Status  delete_edge_points_no_longer_needed( volume, voxel_index,
                                                      voxel_done_flags,
                                                      edge_points )
    volume_struct       *volume;
    voxel_index_struct  *voxel_index;
    unsigned_byte       voxel_done_flags[];
    hash_table_struct   *edge_points;
{
    Status              status;
    int                 axis_index, a1, a2;
    int                 x, y, dx, dy, dz;
    unsigned_byte       get_voxel_done_flag();
    Boolean             all_four_done;
    Boolean             voxel_done[3][3][3];
    voxel_index_struct  indices;
    Status              remove_edge_point();

    status = OK;

    for_inclusive( dx, -1, 1 )
    {
        indices.i[X] = voxel_index->i[X] + dx;
        for_inclusive( dy, -1, 1 )
        {
            indices.i[Y] = voxel_index->i[Y] + dy;
            for_inclusive( dz, -1, 1 )
            {
                indices.i[Z] = voxel_index->i[Z] + dz;

                if( !cube_is_within_volume( volume, indices.i[X],
                             indices.i[Y], indices.i[Z] ) ||
                    get_voxel_done_flag( volume, voxel_done_flags, &indices )
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

                    status = remove_edge_point( volume, edge_points, &indices,
                                                axis_index );
                }
            }
        }
    }

    return( status );
}
