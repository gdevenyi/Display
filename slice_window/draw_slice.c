
#include  <display.h>

typedef  enum  { DIVIDER_INDEX,
                 SLICE1_INDEX,
                 SLICE2_INDEX,
                 SLICE3_INDEX,
                 SLICE4_INDEX,
                 CURSOR1_INDEX1,
                 CURSOR1_INDEX2,
                 CURSOR2_INDEX1,
                 CURSOR2_INDEX2,
                 CURSOR3_INDEX1,
                 CURSOR3_INDEX2,
                 CURSOR4_INDEX1,
                 CURSOR4_INDEX2,
                 TEXT1_INDEX,
                 TEXT2_INDEX,
                 TEXT3_INDEX,
                 TEXT4_INDEX,
                 N_SLICE_MODELS } Slice_model_indices;

typedef enum { X_VOXEL_PROBE_INDEX,
               Y_VOXEL_PROBE_INDEX,
               Z_VOXEL_PROBE_INDEX,
               X_WORLD_PROBE_INDEX,
               Y_WORLD_PROBE_INDEX,
               Z_WORLD_PROBE_INDEX,
               VOXEL_PROBE_INDEX,
               VAL_PROBE_INDEX,
               LABEL_PROBE_INDEX,
               N_READOUT_MODELS     } Slice_readout_indices;

private  void  render_slice_to_pixels(
    display_struct        *slice_window,
    int                   view_index,
    pixels_struct         *pixels );

public  void  initialize_slice_models(
    display_struct    *slice_window )
{
    int            i, view;
    Point          point;
    lines_struct   *lines;
    object_struct  *object;
    model_struct   *model;
    Colour         colour;

    model = get_graphics_model( slice_window, SLICE_MODEL );

    object = create_object( LINES );
    lines = get_lines_ptr( object );

    initialize_lines( lines, Slice_divider_colour );

    fill_Point( point, 0.0, 0.0, 0.0 );
    add_point_to_line( lines, &point );
    add_point_to_line( lines, &point );

    start_new_line( lines );
    add_point_to_line( lines, &point );
    add_point_to_line( lines, &point );

    start_new_line( lines );
    add_point_to_line( lines, &point );
    add_point_to_line( lines, &point );

    add_object_to_model( model, object );

    for_less( view, 0, N_SLICE_VIEWS )
    {
        object = create_object( PIXELS );
        initialize_pixels( get_pixels_ptr(object), 0, 0, 0, 0, 1.0, 1.0,
                           RGB_PIXEL );
        add_object_to_model( model, object );
    }

    for_less( view, 0, N_SLICE_VIEWS )
    {
        /* --- make inner cursor */

        object = create_object( LINES );
        lines = get_lines_ptr( object );
        initialize_lines( lines, Slice_cursor_colour1 );

        lines->n_points = 8;
        lines->n_items = 4;

        ALLOC( lines->points, lines->n_points );
        ALLOC( lines->end_indices, lines->n_items );
        ALLOC( lines->indices, lines->n_points );

        for_less( i, 0, lines->n_items )
            lines->end_indices[i] = 2 * i + 2;

        for_less( i, 0, lines->n_points )
            lines->indices[i] = i;

        add_object_to_model( model, object );

        /* --- make outer cursor */

        object = create_object( LINES );
        lines = get_lines_ptr( object );
        initialize_lines( lines, Slice_cursor_colour2 );

        lines->n_points = 16;
        lines->n_items = 8;

        ALLOC( lines->points, lines->n_points );
        ALLOC( lines->end_indices, lines->n_items );
        ALLOC( lines->indices, lines->n_points );

        for_less( i, 0, lines->n_items )
            lines->end_indices[i] = 2 * i + 2;

        for_less( i, 0, lines->n_points )
            lines->indices[i] = i;

        add_object_to_model( model, object );
    }

    for_less( view, 0, N_SLICE_VIEWS )
    {
        object = create_object( TEXT );

        get_text_ptr(object)->font = (Font_types) Slice_text_font;
        get_text_ptr(object)->size = Slice_text_font_size;
        get_text_ptr(object)->colour = Slice_text_colour;
        add_object_to_model( model, object );
    }

    /* --- initialize readout values */

    model = get_graphics_model( slice_window, SLICE_READOUT_MODEL );

    if( get_model_bitplanes(model) == OVERLAY_PLANES )
        colour = Readout_text_colour;
    else
        colour = Readout_text_rgb_colour;

    for_inclusive( i, 0, N_READOUT_MODELS )
    {
        object = create_object( TEXT );

        get_text_ptr(object)->colour = colour;
        add_object_to_model( model, object );
    }
}

public  void  rebuild_slice_models(
    display_struct    *slice_window )
{
    rebuild_slice_divider( slice_window );
    rebuild_probe( slice_window );
    rebuild_colour_bar( slice_window );

    set_slice_window_all_update( slice_window );
}

public  void  rebuild_slice_divider(
    display_struct    *slice_window )
{
    model_struct   *model;
    Point          *points;
    int            left_panel_width, left_slice_width, right_slice_width;
    int            bottom_slice_height, top_slice_height, text_panel_height;
    int            colour_bar_height;
    int            x_size, y_size;

    model = get_graphics_model(slice_window,SLICE_MODEL);
    points = get_lines_ptr(model->objects[DIVIDER_INDEX])->points;
    G_get_window_size( slice_window->window, &x_size, &y_size );

    get_slice_window_partitions( slice_window,
                                 &left_panel_width, &left_slice_width,
                                 &right_slice_width,
                                 &bottom_slice_height, &top_slice_height,
                                 &text_panel_height, &colour_bar_height );

    fill_Point( points[0], (Real) left_panel_width,               0.0, 0.0 );
    fill_Point( points[1], (Real) left_panel_width, (Real) (y_size-1), 0.0 );

    fill_Point( points[2], (Real) (left_panel_width + left_slice_width),
                                         0.0, 0.0 );
    fill_Point( points[3], (Real) (left_panel_width + left_slice_width),
                           (Real) (y_size-1), 0.0 );

    fill_Point( points[4], (Real) left_panel_width, (Real) bottom_slice_height,
                           0.0 );
    fill_Point( points[5], (Real) (x_size-1),       (Real) bottom_slice_height,
                           0.0 );
}

public  Bitplane_types  get_slice_readout_bitplanes()
{
    if( G_has_overlay_planes() )
        return( (Bitplane_types) Slice_readout_plane );
    else
        return( NORMAL_PLANES );
}

public  void  rebuild_probe(
    display_struct    *slice_window )
{
    model_struct   *model;
    BOOLEAN        active;
    Volume         volume;
    Real           voxel[MAX_DIMENSIONS];
    int            int_voxel[MAX_DIMENSIONS];
    int            label, i, view_index;
    Real           x_world, y_world, z_world;
    text_struct    *text;
    int            sizes[N_DIMENSIONS];
    Real           value, voxel_value;
    int            x_pos, y_pos, x_min, x_max, y_min, y_max;

    active = get_voxel_in_slice_window( slice_window, voxel, &view_index );

    get_text_display_viewport( slice_window, &x_min, &x_max, &y_min, &y_max );

    if( get_slice_window_volume( slice_window, &volume ) )
        get_volume_sizes( volume, sizes );

    convert_voxel_to_world( volume, voxel,
                            &x_world, &y_world, &z_world );

    if( active )
    {
        convert_real_to_int_voxel( N_DIMENSIONS, voxel, int_voxel );

        GET_VOXEL_3D( voxel_value, volume,
                      int_voxel[X], int_voxel[Y], int_voxel[Z] );

        value = CONVERT_VOXEL_TO_VALUE( get_volume(slice_window), voxel_value );

        label = get_volume_label_data( get_label_volume(slice_window),
                                       int_voxel );
    }

    /* --- do slice readout models */

    model = get_graphics_model( slice_window, SLICE_READOUT_MODEL );

    for_less( i, 0, N_READOUT_MODELS )
    {
        x_pos = x_min + Probe_x_pos + (i - X_VOXEL_PROBE_INDEX) * Probe_x_delta;
        y_pos = y_max - Probe_y_pos - (i - X_VOXEL_PROBE_INDEX) * Probe_y_delta
                - (int) ((i - X_VOXEL_PROBE_INDEX) / 3) * Probe_y_pos;

        text = get_text_ptr( model->objects[i] );

        if( active )
        {
            switch( i )
            {
            case X_VOXEL_PROBE_INDEX:
                (void) sprintf( text->string, Slice_probe_x_voxel_format,
                                voxel[X]+1.0 );
                break;
            case Y_VOXEL_PROBE_INDEX:
                (void) sprintf( text->string, Slice_probe_y_voxel_format,
                                voxel[Y]+1.0 );
                break;
            case Z_VOXEL_PROBE_INDEX:
                (void) sprintf( text->string, Slice_probe_z_voxel_format,
                                voxel[Z]+1.0 );
                break;

            case X_WORLD_PROBE_INDEX:
                (void) sprintf( text->string, Slice_probe_x_world_format,
                                x_world );
                break;
            case Y_WORLD_PROBE_INDEX:
                (void) sprintf( text->string, Slice_probe_y_world_format,
                                y_world );
                break;
            case Z_WORLD_PROBE_INDEX:
                (void) sprintf( text->string, Slice_probe_z_world_format,
                                z_world );
                break;
            case VOXEL_PROBE_INDEX:
                (void) sprintf( text->string, Slice_probe_voxel_format,
                                voxel_value );
                break;
            case VAL_PROBE_INDEX:
                (void) sprintf( text->string, Slice_probe_val_format, value );
                break;
            case LABEL_PROBE_INDEX:
                (void) sprintf( text->string, Slice_probe_label_format, label );
                break;
            }
        }
        else
        {
            text->string[0] = (char) 0;
        }

        fill_Point( text->origin, x_pos, y_pos, 0.0 );
    }

    set_update_required( slice_window, get_slice_readout_bitplanes() );
}

private  void  get_cursor_size(
    int    slice_index,
    Real   *hor_start,
    Real   *hor_end,
    Real   *vert_start,
    Real   *vert_end )
{
    switch( slice_index )
    {
    case 0:
        *hor_start = Cursor_hor_start_0;
        *hor_end = Cursor_hor_end_0;
        *vert_start = Cursor_vert_start_0;
        *vert_end = Cursor_vert_end_0;
        break;

    case 1:
        *hor_start = Cursor_hor_start_1;
        *hor_end = Cursor_hor_end_1;
        *vert_start = Cursor_vert_start_1;
        *vert_end = Cursor_vert_end_1;
        break;

    case 2:
        *hor_start = Cursor_hor_start_2;
        *hor_end = Cursor_hor_end_2;
        *vert_start = Cursor_vert_start_2;
        *vert_end = Cursor_vert_end_2;
        break;

    case 3:
        *hor_start = Cursor_hor_start_3;
        *hor_end = Cursor_hor_end_3;
        *vert_start = Cursor_vert_start_3;
        *vert_end = Cursor_vert_end_3;
        break;
    }
}

private  void  rebuild_cursor(
    display_struct    *slice_window,
    int               view_index )
{
    model_struct   *model;
    int            c, x_index, y_index, axis;
    Real           x_left, x_right, y_bottom, y_top, dx, dy;
    Real           x_centre, y_centre;
    Real           tmp_voxel[N_DIMENSIONS];
    lines_struct   *lines1, *lines2;
    Real           current_voxel[N_DIMENSIONS];
    int            x_min, x_max, y_min, y_max;
    Real           hor_pixel_start, hor_pixel_end;
    Real           vert_pixel_start, vert_pixel_end;

    model = get_graphics_model(slice_window,SLICE_MODEL);

    lines1 = get_lines_ptr( model->objects[CURSOR1_INDEX1+2*view_index] );
    lines2 = get_lines_ptr( model->objects[CURSOR1_INDEX1+2*view_index+1] );

    get_current_voxel( slice_window, current_voxel );

    if( slice_has_ortho_axes( slice_window, view_index,
                              &x_index, &y_index, &axis ) )
    {
        for_less( c, 0, N_DIMENSIONS )
            tmp_voxel[c] = ROUND( current_voxel[c] );

        current_voxel[x_index] += 0.5;
        convert_voxel_to_pixel( slice_window, view_index, current_voxel,
                                &x_right, &y_centre );
        current_voxel[x_index] -= 0.5;

        current_voxel[x_index] -= 0.5;
        convert_voxel_to_pixel( slice_window, view_index, current_voxel,
                                &x_left, &y_centre );
        current_voxel[x_index] += 0.5;

        current_voxel[y_index] += 0.5;
        convert_voxel_to_pixel( slice_window, view_index, current_voxel,
                                &x_centre, &y_top );
        current_voxel[y_index] -= 0.5;

        current_voxel[y_index] -= 0.5;
        convert_voxel_to_pixel( slice_window, view_index, current_voxel,
                                &x_centre, &y_bottom );
        current_voxel[y_index] += 0.5;
    }
    else
    {
        convert_voxel_to_pixel( slice_window, view_index, current_voxel,
                                &x_centre, &y_centre );
        x_left = x_centre;
        x_right = x_centre;
        y_bottom = y_centre;
        y_top = y_centre;
    }

    get_cursor_size( view_index, &hor_pixel_start, &hor_pixel_end,
                     &vert_pixel_start, &vert_pixel_end );

    get_slice_viewport( slice_window, view_index,
                        &x_min, &x_max, &y_min, &y_max );

    if( x_centre < x_min )
    {
        dx = x_min - x_centre;
        x_centre = x_min;
        x_left += dx;
        x_right += dx;
    }
    else if( x_centre > x_max )
    {
        dx = x_max - x_centre;
        x_centre = x_max;
        x_left += dx;
        x_right += dx;
    }

    if( y_centre < y_min )
    {
        dy = y_min - y_centre;
        y_centre = y_min;
        y_top += dy;
        y_bottom += dy;
    }
    else if( y_centre > y_max )
    {
        dy = y_max - y_centre;
        y_centre = y_max;
        y_top += dy;
        y_bottom += dy;
    }

    fill_Point( lines1->points[0], x_right + hor_pixel_start, y_centre, 0.0 );
    fill_Point( lines1->points[1], x_right + hor_pixel_end, y_centre, 0.0 );
    fill_Point( lines1->points[2], x_left - hor_pixel_start, y_centre, 0.0 );
    fill_Point( lines1->points[3], x_left - hor_pixel_end, y_centre, 0.0 );
    fill_Point( lines1->points[4], x_centre, y_top + vert_pixel_start, 0.0 );
    fill_Point( lines1->points[5], x_centre, y_top + vert_pixel_end, 0.0 );
    fill_Point( lines1->points[6], x_centre, y_bottom - vert_pixel_start, 0.0 );
    fill_Point( lines1->points[7], x_centre, y_bottom - vert_pixel_end, 0.0 );

    fill_Point( lines2->points[0], x_right + hor_pixel_start, y_centre-1.0,0.0);
    fill_Point( lines2->points[1], x_right + hor_pixel_end, y_centre-1.0, 0.0 );
    fill_Point( lines2->points[2], x_right + hor_pixel_start, y_centre+1.0,0.0);
    fill_Point( lines2->points[3], x_right + hor_pixel_end, y_centre+1.0, 0.0 );

    fill_Point( lines2->points[4], x_left - hor_pixel_start, y_centre-1.0, 0.0);
    fill_Point( lines2->points[5], x_left - hor_pixel_end, y_centre-1.0, 0.0 );
    fill_Point( lines2->points[6], x_left - hor_pixel_start, y_centre+1.0, 0.0);
    fill_Point( lines2->points[7], x_left - hor_pixel_end, y_centre+1.0, 0.0 );

    fill_Point( lines2->points[8], x_centre-1.0, y_top + vert_pixel_start, 0.0);
    fill_Point( lines2->points[9], x_centre-1.0, y_top + vert_pixel_end, 0.0 );
    fill_Point( lines2->points[10],x_centre+1.0, y_top + vert_pixel_start, 0.0);
    fill_Point( lines2->points[11],x_centre+1.0, y_top + vert_pixel_end, 0.0 );

    fill_Point( lines2->points[12],x_centre-1.0, y_bottom-vert_pixel_start,0.0);
    fill_Point( lines2->points[13],x_centre-1.0, y_bottom-vert_pixel_end,0.0);
    fill_Point( lines2->points[14],x_centre+1.0, y_bottom-vert_pixel_start,0.0);
    fill_Point( lines2->points[15],x_centre+1.0, y_bottom-vert_pixel_end,0.0);
}

public  void  rebuild_cursors(
    display_struct    *slice_window )
{
    int   view;

    for_less( view, 0, N_SLICE_VIEWS )
        rebuild_cursor( slice_window, view );
}

public  object_struct  *get_slice_pixels_object(
    display_struct    *slice_window,
    int               view_index )
{
    model_struct   *model;

    model = get_graphics_model(slice_window,SLICE_MODEL);

    return( model->objects[SLICE1_INDEX+view_index] );
}

public  void  rebuild_slice_pixels(
    display_struct    *slice_window,
    int               view_index )
{
    BOOLEAN        visibility;
    model_struct   *model;
    object_struct  *pixels_object;
    pixels_struct  *pixels;
    int            axis_index, x_index, y_index;
    int            x_min, x_max, y_min, y_max;
    text_struct    *text;
    char           *format;
    int            x_pos, y_pos;
    Real           current_voxel[N_DIMENSIONS];

    pixels_object = get_slice_pixels_object(slice_window,view_index);
    pixels = get_pixels_ptr( pixels_object );

    visibility = get_slice_visibility( slice_window, view_index );

    set_object_visibility( pixels_object, visibility );

    if( visibility )
        render_slice_to_pixels( slice_window, view_index, pixels );

    model = get_graphics_model(slice_window,SLICE_MODEL);

    if( slice_has_ortho_axes( slice_window, view_index,
                              &x_index, &y_index, &axis_index ) )
    {
        set_object_visibility( model->objects[TEXT1_INDEX+view_index], TRUE );

        text = get_text_ptr( model->objects[TEXT1_INDEX+view_index] );

        switch( axis_index )
        {
        case X:  format = Slice_index_x_format;  break;
        case Y:  format = Slice_index_y_format;  break;
        case Z:  format = Slice_index_z_format;  break;
        }

        get_current_voxel( slice_window, current_voxel );

        (void) sprintf( text->string, format, current_voxel[axis_index] + 1.0 );

        get_slice_viewport( slice_window,
                            view_index, &x_min, &x_max, &y_min, &y_max );

        x_pos = x_min + (int) Point_x(Slice_index_offset);
        y_pos = y_min + (int) Point_y(Slice_index_offset);

        fill_Point( text->origin, x_pos, y_pos, 0.0 );
    }
    else
        set_object_visibility( model->objects[TEXT1_INDEX+view_index], FALSE );

    rebuild_cursors( slice_window );
}

#define  MAX_LABELS   256

private  void  render_slice_to_pixels(
    display_struct        *slice_window,
    int                   view_index,
    pixels_struct         *pixels )
{
    Volume                volume, label_volume;
    int                   x_size, y_size;
    int                   n_alloced;
    Real                  x_trans, y_trans, x_scale, y_scale;
    Real                  origin[MAX_DIMENSIONS];
    Real                  x_axis[MAX_DIMENSIONS], y_axis[MAX_DIMENSIONS];
    int                   x_min, x_max, y_min, y_max;

    if( pixels->x_size > 0 && pixels->y_size > 0 )
        delete_pixels( pixels );

    n_alloced = 0;

    x_trans = slice_window->slice.slice_views[view_index].x_trans;
    y_trans = slice_window->slice.slice_views[view_index].y_trans;
    x_scale = slice_window->slice.slice_views[view_index].x_scaling;
    y_scale = slice_window->slice.slice_views[view_index].y_scaling;

    get_slice_viewport( slice_window, view_index,
                        &x_min, &x_max, &y_min, &y_max );

    volume = get_volume( slice_window );
    label_volume = get_label_volume( slice_window );

    get_slice_plane( slice_window, view_index, origin, x_axis, y_axis );

    if( slice_window->slice.display_labels &&
        label_volume->data != (void *) NULL )
    {
        create_volume_slice(

                    label_volume, NEAREST_NEIGHBOUR, 0.0,
                    origin, x_axis, y_axis,
                    x_trans, y_trans, x_scale, y_scale,

                    volume,
                    slice_window->slice.slice_views[view_index].filter_type,
                    slice_window->slice.slice_views[view_index].filter_width,
                    origin, x_axis, y_axis,
                    x_trans, y_trans, x_scale, y_scale,

                    x_max - x_min + 1, y_max - y_min + 1,
                    RGB_PIXEL, FALSE, (unsigned short **) NULL,
                    slice_window->slice.colour_tables,
                    make_rgba_Colour( 0, 0, 0, 0 ),
                    &n_alloced, pixels );
    }
    else
    {
        create_volume_slice(
                    volume,
                    slice_window->slice.slice_views[view_index].filter_type,
                    slice_window->slice.slice_views[view_index].filter_width,
                    origin, x_axis, y_axis,
                    x_trans, y_trans, x_scale, y_scale,
                    (Volume) NULL, NEAREST_NEIGHBOUR, 0.0,
                    (Real *) 0, (Real *) 0, (Real *) 0,
                    0.0, 0.0, 0.0, 0.0,
                    x_max - x_min + 1, y_max - y_min + 1,
                    RGB_PIXEL, FALSE, (unsigned short **) NULL,
                    slice_window->slice.colour_tables,
                    make_rgba_Colour( 0, 0, 0, 0 ),
                    &n_alloced, pixels );

    }

    pixels->x_position += x_min;
    pixels->y_position += y_min;

    x_size = pixels->x_size;
    y_size = pixels->y_size;

    /* --- now blend in the talaiarch atlas */

    if( x_size > 0 && y_size > 0 )
    {
        Real  v1[N_DIMENSIONS], v2[N_DIMENSIONS];
        Real  dx, dy;
        int   sizes[N_DIMENSIONS];
        int   x_index, y_index, axis_index;

        if( !slice_has_ortho_axes( slice_window, view_index,
                                   &x_index, &y_index, &axis_index ) )
            return;

        (void) convert_slice_pixel_to_voxel( volume,
                        pixels->x_position - x_min, pixels->y_position - y_min,
                        origin, x_axis, y_axis,
                        x_trans, y_trans, x_scale, y_scale, v1 );
        (void) convert_slice_pixel_to_voxel( volume,
                        pixels->x_position+1 - x_min, pixels->y_position-y_min,
                        origin, x_axis, y_axis,
                        x_trans, y_trans, x_scale, y_scale, v2 );

        dx = v2[x_index] - v1[x_index];

        (void) convert_slice_pixel_to_voxel( volume,
                        pixels->x_position - x_min, pixels->y_position+1-y_min,
                        origin, x_axis, y_axis,
                        x_trans, y_trans, x_scale, y_scale, v2 );

        dy = v2[y_index] - v1[y_index];

        get_volume_sizes( volume, sizes );

        blend_in_atlas( &slice_window->slice.atlas,
                        pixels->data.pixels_rgb,
                        x_size, y_size,
                        v1, x_index, y_index, axis_index,
                        dx, dy,
                        sizes[x_index], sizes[y_index] );
    }
}
