
#include  <def_graphics.h>
#include  <def_globals.h>

public  Status  initialize_cursor( graphics )
    graphics_struct   *graphics;
{
    Status          status;
    Status          reset_cursor();
    Status          initialize_cursor_plane_outline();

    graphics->three_d.cursor.box_size[X] = 1.0;
    graphics->three_d.cursor.box_size[Y] = 1.0;
    graphics->three_d.cursor.box_size[Z] = 1.0;
    graphics->three_d.cursor.axis_size = Cursor_axis_size;

    status = reset_cursor( graphics );

    if( status == OK )
        status = initialize_cursor_plane_outline( graphics );

    return( status );
}

public  Status  reset_cursor( graphics )
    graphics_struct   *graphics;
{
    Status          status;
    Status          rebuild_cursor_icon();
    Real            size_of_domain();
    void            update_cursor();

    graphics->three_d.cursor.origin = graphics->three_d.centre_of_objects;

    status = rebuild_cursor_icon( graphics );

    graphics->models[CURSOR_MODEL]->visibility = ON;

    update_cursor( graphics );

    return( status );
}

public  Status  update_cursor_size( graphics )
    graphics_struct   *graphics;
{
    Status          rebuild_cursor_icon();

    return( rebuild_cursor_icon( graphics ) );
}

public  void  update_cursor( graphics )
    graphics_struct   *graphics;
{
    void           make_origin_transform();
    model_struct   *model;
    model_struct   *get_graphics_model();

    model = get_graphics_model( graphics, CURSOR_MODEL );

    make_origin_transform( &graphics->three_d.cursor.origin,
                           &model->transform );

    ++graphics->models_changed_id;
}
