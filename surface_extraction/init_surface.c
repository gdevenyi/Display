
#include  <def_graphics.h>
#include  <def_globals.h>
#include  <def_marching_cubes.h>
#include  <def_splines.h>
#include  <def_bitlist.h>
#include  <def_files.h>

static    Status   clear_surface_extraction();

public  Status  initialize_surface_extraction( graphics )
    graphics_struct    *graphics;
{
    Status                      status;
    surface_extraction_struct   *surface_extraction;
    Status                      create_object();
    Status                      add_object_to_model();
    object_struct               *object;
    void                        install_surface_extraction();

    surface_extraction = &graphics->three_d.surface_extraction;

    status = create_object( &object, POLYGONS );

    if( status == OK )
    {
        status = add_object_to_model(
                     graphics->models[THREED_MODEL]->ptr.model, object );
    }

    if( status == OK )
    {
        surface_extraction->polygons = object->ptr.polygons;

        surface_extraction->polygons->colour_flag = ONE_COLOUR;

        ALLOC( status, surface_extraction->polygons->colours, 1 );
    }

    if( status == OK )
        status = create_bitlist( 0, &surface_extraction->voxels_queued );

    if( status == OK )
    {
        install_surface_extraction( graphics );

        status = clear_surface_extraction( graphics );
    }

    return( status );
}

private  Status  clear_surface_extraction( graphics )
    graphics_struct    *graphics;
{
    Status                      status;
    void                        empty_polygons_struct();
    Status                      initialize_edge_points();
    void                        initialize_voxel_queue();
    void                        clear_voxel_flags();
    void                        clear_voxel_done_flags();
    surface_extraction_struct   *surface_extraction;
    volume_struct               *volume;
    Status                      delete_polygons();

    surface_extraction = &graphics->three_d.surface_extraction;

    surface_extraction->extraction_in_progress = FALSE;
    surface_extraction->isovalue_selected = FALSE;
    surface_extraction->n_voxels_with_surface = 0;

    status = initialize_edge_points( &surface_extraction->edge_points );

    if( status == OK )
    {
        status = delete_polygons( surface_extraction->polygons );
    }

    if( status == OK )
    {
        initialize_voxel_queue( &surface_extraction->voxels_to_do );

        empty_polygons_struct( surface_extraction->polygons );

        ALLOC( status, surface_extraction->polygons->colours, 1 );

        surface_extraction->polygons->colours[0] = Extracted_surface_colour;
        surface_extraction->polygons->surfprop = Default_surface_property;

        clear_voxel_flags( &surface_extraction->voxels_queued );

        if( get_slice_window_volume( graphics, &volume ) )
        {
            clear_voxel_done_flags( surface_extraction->voxel_done_flags,
                                    get_n_voxels(volume) );
        }
        else
        {
            surface_extraction->voxel_done_flags = (unsigned_byte *) 0;
        }
    }

    return( status );
}

private  Status  free_surface_extraction( graphics )
    graphics_struct   *graphics;
{
    Status                      status;
    Status                      delete_edge_points();
    Status                      delete_voxel_queue();
    surface_extraction_struct   *surface_extraction;

    surface_extraction = &graphics->three_d.surface_extraction;

    status = delete_edge_points( &surface_extraction->edge_points );

    if( status == OK )
    {
        status = delete_voxel_queue( &surface_extraction->voxels_to_do );
    }

    return( status );
}

public  Status  delete_surface_extraction( graphics )
    graphics_struct   *graphics;
{
    Status                      status;
    Status                      delete_voxel_flags();
    Status                      delete_voxel_done_flags();
    surface_extraction_struct   *surface_extraction;

    surface_extraction = &graphics->three_d.surface_extraction;

    status = free_surface_extraction( graphics );

    if( status == OK )
    {
        status = delete_voxel_done_flags( surface_extraction->voxel_done_flags);
    }

    if( status == OK )
    {
        status = delete_voxel_flags( &surface_extraction->voxels_queued );
    }

    return( status );
}

public  Status  reset_surface_extraction( graphics )
    graphics_struct   *graphics;
{
    Status    status;
    Status    free_surface_extraction();

    status = free_surface_extraction( graphics );

    if( status == OK )
    {
        status = clear_surface_extraction( graphics );
    }

    return( status );
}

public  void  set_isosurface_value( surface_extraction )
    surface_extraction_struct   *surface_extraction;
{
    Real   value;

    if( !surface_extraction->extraction_in_progress )
    {
        PRINT( "Enter isosurface value: " );
        if( input_real( stdin, &value ) == OK && value >= 0.0 )
        {
            if( value == (Real) ((int) value) )
            {
                value = value + 0.0001;
            }

            surface_extraction->isovalue = value;
            surface_extraction->isovalue_selected = TRUE;
        }

        (void) input_newline( stdin );
    }
}

public  void  check_if_isosurface_value_set( surface_extraction )
    surface_extraction_struct   *surface_extraction;
{
    void    set_isosurface_value();

    if( !surface_extraction->isovalue_selected )
    {
        set_isosurface_value( surface_extraction );
    }
}

public  Boolean  get_isosurface_value( graphics, value )
    graphics_struct    *graphics;
    Real               *value;
{
    if( graphics->three_d.surface_extraction.isovalue_selected )
    {
        *value = graphics->three_d.surface_extraction.isovalue;
    }

    return( graphics->three_d.surface_extraction.isovalue_selected );
}

public  void  start_surface_extraction( graphics )
    graphics_struct    *graphics;
{
    surface_extraction_struct   *surface_extraction;

    surface_extraction = &graphics->three_d.surface_extraction;

    if( !surface_extraction->extraction_in_progress &&
        surface_extraction->isovalue_selected )
    {
        surface_extraction->extraction_in_progress = TRUE;
    }
}

public  void  stop_surface_extraction( graphics )
    graphics_struct    *graphics;
{
    if( graphics->three_d.surface_extraction.extraction_in_progress )
    {
        graphics->three_d.surface_extraction.extraction_in_progress = FALSE;
    }
}

public  int  get_n_voxels( volume )
    volume_struct  *volume;
{
    int   n_voxels;
    int   nx, ny, nz;
    void  get_volume_size();

    if( volume != (volume_struct *) 0 )
    {
        get_volume_size( volume, &nx, &ny, &nz );

        n_voxels = (nx - 1) * (ny - 1) * (nz - 1);
    }
    else
    {
        n_voxels = 0;
    }

    return( n_voxels );
}
