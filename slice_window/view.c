
#include  <display.h>

public  void  initialize_slice_window_view(
    display_struct    *slice_window )
{
    Volume           volume;
    int              axis, view;
    int              size[N_DIMENSIONS];
    Real             thickness[N_DIMENSIONS];

    slice_window->slice.slice_views[0].axis_map[0]  = Slice_view1_axis1;
    slice_window->slice.slice_views[0].axis_map[1]  = Slice_view1_axis2;
    slice_window->slice.slice_views[0].axis_map[2]  = Slice_view1_axis3;

    slice_window->slice.slice_views[1].axis_map[0]  = Slice_view2_axis1;
    slice_window->slice.slice_views[1].axis_map[1]  = Slice_view2_axis2;
    slice_window->slice.slice_views[1].axis_map[2]  = Slice_view2_axis3;

    slice_window->slice.slice_views[2].axis_map[0]  = Slice_view3_axis1;
    slice_window->slice.slice_views[2].axis_map[1]  = Slice_view3_axis2;
    slice_window->slice.slice_views[2].axis_map[2]  = Slice_view3_axis3;

    if( get_slice_window_volume( slice_window, &volume ) )
    {
        get_volume_sizes( volume, size );
        get_volume_separations( volume, thickness );

        for_less( view, 0, N_DIMENSIONS )
        {
            axis = slice_window->slice.slice_views[view].axis_map[Z];
            slice_window->slice.slice_index[view] = size[axis] / 2.0;
        }

        for_less( view, 0, N_DIMENSIONS )
            reset_slice_view( slice_window, view );

        slice_window->associated[THREE_D_WINDOW]->three_d.cursor.box_size[X] =
                          ABS( thickness[X] );
        slice_window->associated[THREE_D_WINDOW]->three_d.cursor.box_size[Y] =
                          ABS( thickness[Y] );
        slice_window->associated[THREE_D_WINDOW]->three_d.cursor.box_size[Z] =
                          ABS( thickness[Z] );

        update_cursor_size( slice_window->associated[THREE_D_WINDOW] );
    }
}

public  void  reset_slice_view(
    display_struct    *slice_window,
    int               view )
{
    Volume  volume;
    int     x_min, x_max, y_min, y_max;
    Real    origin[MAX_DIMENSIONS];
    Real    x_axis[MAX_DIMENSIONS], y_axis[MAX_DIMENSIONS];;

    if( get_slice_window_volume( slice_window, &volume ) )
    {
        get_slice_viewport( slice_window, view,
                            &x_min, &x_max, &y_min, &y_max );

        get_slice_plane( slice_window, view, origin, x_axis, y_axis );

        fit_volume_slice_to_viewport( volume, origin, x_axis, y_axis,
                x_max - x_min + 1, y_max - y_min + 1,
                Slice_fit_oversize,
                &slice_window->slice.slice_views[view].x_trans,
                &slice_window->slice.slice_views[view].y_trans,
                &slice_window->slice.slice_views[view].x_scaling,
                &slice_window->slice.slice_views[view].y_scaling,
                &slice_window->slice.slice_views[view].used_viewport_x_size,
                &slice_window->slice.slice_views[view].used_viewport_y_size );

         slice_window->slice.slice_views[view].prev_viewport_x_size =
                                                          (x_max - x_min + 1);
         slice_window->slice.slice_views[view].prev_viewport_y_size =
                                                          (y_max - y_min + 1);
    }
}

public  void  resize_slice_view(
    display_struct    *slice_window,
    int               view )
{
    int            x_min, x_max, y_min, y_max;
    Volume         volume;

    if( get_slice_window_volume( slice_window, &volume ) )
    {
        get_slice_viewport( slice_window, view,
                            &x_min, &x_max, &y_min, &y_max );

        resize_volume_slice(
             slice_window->slice.slice_views[view].prev_viewport_x_size,
             slice_window->slice.slice_views[view].prev_viewport_y_size,
             slice_window->slice.slice_views[view].used_viewport_x_size,
             slice_window->slice.slice_views[view].used_viewport_y_size,
             x_max - x_min + 1, y_max - y_min + 1,
             &slice_window->slice.slice_views[view].x_trans,
             &slice_window->slice.slice_views[view].y_trans,
             &slice_window->slice.slice_views[view].x_scaling,
             &slice_window->slice.slice_views[view].y_scaling,
             &slice_window->slice.slice_views[view].used_viewport_x_size,
             &slice_window->slice.slice_views[view].used_viewport_y_size );

         slice_window->slice.slice_views[view].prev_viewport_x_size =
                                                          (x_max - x_min + 1);
         slice_window->slice.slice_views[view].prev_viewport_y_size =
                                                          (y_max - y_min + 1);
    }
}

public  void  scale_slice_view(
    display_struct    *slice_window,
    int               view,
    Real              scale_factor )
{
    Volume volume;
    int    x_min, x_max, y_min, y_max;

    if( get_slice_window_volume( slice_window, &volume ) )
    {
        get_slice_viewport( slice_window, view,
                            &x_min, &x_max, &y_min, &y_max );

        scale_slice_about_viewport_centre( scale_factor,
             x_max - x_min + 1, y_max - y_min + 1,
             &slice_window->slice.slice_views[view].x_trans,
             &slice_window->slice.slice_views[view].y_trans,
             &slice_window->slice.slice_views[view].x_scaling,
             &slice_window->slice.slice_views[view].y_scaling );
    }
}

public  void  translate_slice_view(
    display_struct    *slice_window,
    int               view,
    Real              dx,
    Real              dy )
{
    Volume volume;

    if( get_slice_window_volume( slice_window, &volume ) )
    {
        slice_window->slice.slice_views[view].x_trans += dx;
        slice_window->slice.slice_views[view].y_trans += dy;
    }
}

public  BOOLEAN  find_slice_view_mouse_is_in(
    display_struct    *display,
    int               x_pixel,
    int               y_pixel,
    int               *view_index )
{
    BOOLEAN          found;
    int              c;
    int              x_min, x_max, y_min, y_max;
    display_struct   *slice_window;

    found = FALSE;

    if( get_slice_window(display,&slice_window) )
    {
        for_less( c, 0, N_DIMENSIONS )
        {
            get_slice_viewport( slice_window, c,
                                &x_min, &x_max, &y_min, &y_max );

            if( x_pixel >= x_min && x_pixel <= x_max &&
                y_pixel >= y_min && y_pixel <= y_max )
            {
                *view_index = c;
                found = TRUE;

                break;
            }
        }
    }

    return( found );
}

public  BOOLEAN  convert_pixel_to_voxel(
    display_struct    *display,
    int               x_pixel,
    int               y_pixel,
    Real              voxel[],
    int               *view_index )
{
    BOOLEAN           found;
    Volume            volume;
    display_struct    *slice_window;
    int               x_min, x_max, y_min, y_max;
    Real              origin[MAX_DIMENSIONS];
    Real              x_axis[MAX_DIMENSIONS], y_axis[MAX_DIMENSIONS];;

    found = FALSE;

    if( find_slice_view_mouse_is_in( display, x_pixel, y_pixel, view_index ) &&
        get_slice_window_volume( display, &volume ) &&
        get_slice_window( display, &slice_window ) )
    {
        get_slice_viewport( slice_window, *view_index,
                            &x_min, &x_max, &y_min, &y_max );
        get_slice_plane( slice_window, *view_index, origin, x_axis, y_axis );

        x_pixel -= x_min;
        y_pixel -= y_min;

        found = convert_slice_pixel_to_voxel( volume, x_pixel, y_pixel,
                    origin, x_axis, y_axis,
                    slice_window->slice.slice_views[*view_index].x_trans,
                    slice_window->slice.slice_views[*view_index].y_trans,
                    slice_window->slice.slice_views[*view_index].x_scaling,
                    slice_window->slice.slice_views[*view_index].y_scaling,
                    voxel );
    }

    return( found );
}

public  void  convert_voxel_to_pixel(
    display_struct    *display,
    int               view_index,
    Real              voxel[],
    int               *x_pixel,
    int               *y_pixel )
{
    Volume            volume;
    display_struct    *slice_window;
    int               x_min, x_max, y_min, y_max;
    Real              x_real_pixel, y_real_pixel;
    Real              origin[MAX_DIMENSIONS];
    Real              x_axis[MAX_DIMENSIONS], y_axis[MAX_DIMENSIONS];;

    if( get_slice_window( display, &slice_window ) &&
        get_slice_window_volume( display, &volume ) )
    {
        get_slice_plane( slice_window, view_index, origin, x_axis, y_axis );

        convert_voxel_to_slice_pixel( volume, voxel, origin, x_axis, y_axis,
                        slice_window->slice.slice_views[view_index].x_trans,
                        slice_window->slice.slice_views[view_index].y_trans,
                        slice_window->slice.slice_views[view_index].x_scaling,
                        slice_window->slice.slice_views[view_index].y_scaling,
                        &x_real_pixel, &y_real_pixel );

        get_slice_viewport( display, view_index,
                            &x_min, &x_max, &y_min, &y_max );

        *x_pixel = ROUND( x_real_pixel ) + x_min;
        *y_pixel = ROUND( y_real_pixel ) + y_min;
    }
    else
    {
        HANDLE_INTERNAL_ERROR( "convert_voxel_to_pixel" );
    }
}

public  BOOLEAN  get_voxel_corresponding_to_point(
    display_struct    *display,
    Point             *point,
    Real              voxel[] )
{
    Volume          volume;
    BOOLEAN         converted;

    converted = FALSE;

    if( get_slice_window_volume( display, &volume ) )
    {
        convert_world_to_voxel( volume,
                            Point_x(*point), Point_y(*point), Point_z(*point),
                            voxel );

        converted = voxel_is_within_volume( volume, voxel );
    }

    return( converted );
}

public  void  get_slice_viewport(
    display_struct    *display,
    int               view_index,
    int               *x_min,
    int               *x_max,
    int               *y_min,
    int               *y_max )
{
    int  x_size, y_size;

    G_get_window_size( display->window, &x_size, &y_size );

    switch( view_index )
    {
    case 0:
        *x_min = Slice_divider_left;
        *x_max = display->slice.x_split-1-Slice_divider_right;
        *y_min = display->slice.y_split+1+Slice_divider_bottom;
        *y_max = y_size-Slice_divider_top;
        break;

    case 1:
        *x_min = display->slice.x_split+1+Slice_divider_left;
        *x_max = x_size-Slice_divider_right;
        *y_min = display->slice.y_split+1+Slice_divider_bottom;
        *y_max = y_size-Slice_divider_top;
        break;

    case 2:
        *x_min = Slice_divider_left;
        *x_max = display->slice.x_split-1-Slice_divider_right;
        *y_min = Slice_divider_bottom;
        *y_max = display->slice.y_split-1-Slice_divider_top;
        break;

    default:
        *x_min = display->slice.x_split+1+Slice_divider_left;
        *x_max = x_size-Slice_divider_right;
        *y_min = Slice_divider_bottom;
        *y_max = display->slice.y_split-1-Slice_divider_top;
    }
}

public  BOOLEAN  get_voxel_in_slice_window(
    display_struct    *display,
    Real              voxel[],
    int               *view_index )
{
    display_struct    *slice_window;
    int               x_mouse, y_mouse;
    BOOLEAN           found;

    slice_window = display->associated[SLICE_WINDOW];

    (void) G_get_mouse_position( slice_window->window, &x_mouse, &y_mouse );

    found = convert_pixel_to_voxel( slice_window, x_mouse, y_mouse, voxel,
                                    view_index );

    return( found );
}

public  BOOLEAN  get_voxel_in_three_d_window(
    display_struct    *display,
    Real              voxel[] )
{
    BOOLEAN          found;
    object_struct    *object;
    int              object_index;
    Point            intersection_point;
    display_struct   *slice_window;

    found = FALSE;

    if( get_mouse_scene_intersection( display, &object, &object_index,
                                      &intersection_point ) &&
        get_slice_window( display, &slice_window ) )
    {
        found = get_voxel_corresponding_to_point( slice_window,
                                                  &intersection_point,
                                                  voxel );
    }

    return( found );
}

public  BOOLEAN  get_voxel_under_mouse(
    display_struct    *display,
    Real              voxel[],
    int               *view_index )
{
    display_struct    *three_d, *slice_window;
    BOOLEAN           found;

    three_d = display->associated[THREE_D_WINDOW];

    if( !get_slice_window(three_d,&slice_window) )
        found = FALSE;
    else if( G_is_mouse_in_window( slice_window->window ) )
    {
        found = get_voxel_in_slice_window( display, voxel, view_index );
    }
    else if( G_is_mouse_in_window( three_d->window ) )
    {
        found = get_voxel_in_three_d_window( three_d, voxel );
        *view_index = 2;
    }
    else
    {
        found = FALSE;
    }

    return( found );
}

public  void  get_current_voxel(
    display_struct    *slice_window,
    Real              voxel[] )
{
    voxel[X] = slice_window->slice.slice_index[X];
    voxel[Y] = slice_window->slice.slice_index[Y];
    voxel[Z] = slice_window->slice.slice_index[Z];
}

public  BOOLEAN  set_current_voxel(
    display_struct    *slice_window,
    Real              voxel[] )
{
    BOOLEAN           changed;
    int               i, view, sizes[MAX_DIMENSIONS];
    int               axis, x_index, y_index;

    changed = FALSE;

    get_volume_sizes( get_volume(slice_window), sizes );

    for_less( i, 0, N_DIMENSIONS )
    {
        if( voxel[i] != slice_window->slice.slice_index[i] )
        {
            slice_window->slice.slice_index[i] = voxel[i];

            for_less( view, 0, N_DIMENSIONS )
            {
                if( !slice_has_ortho_axes( slice_window, view,
                                           &x_index, &y_index, &axis ) ||
                    axis == i )
                {
                    set_slice_window_update( slice_window, view );
                }
            }

            changed = TRUE;
        }
    }

    return( changed );
}

public  void  get_slice_axes(
    display_struct   *slice_window,
    int              view_index,
    int              *x_index,
    int              *y_index,
    int              *axis_index )
{
    *x_index = slice_window->slice.slice_views[view_index].axis_map[X];
    *y_index = slice_window->slice.slice_views[view_index].axis_map[Y];
    *axis_index = slice_window->slice.slice_views[view_index].axis_map[Z];
}

public  void  get_slice_plane(
    display_struct   *slice_window,
    int              view_index,
    Real             origin[],
    Real             x_axis[],
    Real             y_axis[] )
{
    int    x_index, y_index, axis;
    Real   voxel_indices[MAX_DIMENSIONS];

    get_slice_axes( slice_window, view_index, &x_index, &y_index, &axis );
    get_current_voxel( slice_window, voxel_indices );

    origin[X] = 0.0;
    origin[Y] = 0.0;
    origin[Z] = 0.0;
    origin[axis] = voxel_indices[axis];

    x_axis[X] = 0.0;
    x_axis[Y] = 0.0;
    x_axis[Z] = 0.0;
    x_axis[x_index] = 1.0;

    y_axis[X] = 0.0;
    y_axis[Y] = 0.0;
    y_axis[Z] = 0.0;
    y_axis[y_index] = 1.0;
}

public  void  get_slice_perp_axis(
    display_struct   *slice_window,
    int              view_index,
    Real             perp_axis[N_DIMENSIONS] )
{
    Real     origin[N_DIMENSIONS];
    Real     x_axis[N_DIMENSIONS];
    Real     y_axis[N_DIMENSIONS];
    int      c, a1, a2;

    get_slice_plane( slice_window, view_index, origin, x_axis, y_axis );

    for_less( c, 0, N_DIMENSIONS )
    {
        a1 = (c + 1) % N_DIMENSIONS;
        a2 = (c + 2) % N_DIMENSIONS;
        perp_axis[c] = x_axis[a1] * y_axis[a2] - x_axis[a2] * y_axis[a1];
    }
}

public  BOOLEAN  get_axis_index_under_mouse(
    display_struct   *display,
    int              *axis_index )
{
    BOOLEAN          found;
    int              view_index;
    display_struct   *slice_window;

    found = get_slice_view_index_under_mouse( display, &view_index );

    if( found )
    {
        slice_window = display->associated[SLICE_WINDOW];

        *axis_index =
             slice_window->slice.slice_views[view_index].axis_map[Z];
    }

    return( found );
}

public  BOOLEAN  slice_has_ortho_axes(
    display_struct   *slice_window,
    int              view_index,
    int              *x_index,
    int              *y_index,
    int              *axis_index )
{
    Real     origin[N_DIMENSIONS];
    Real     x_axis[N_DIMENSIONS];
    Real     y_axis[N_DIMENSIONS];
    int      c;

    get_slice_plane( slice_window, view_index, origin, x_axis, y_axis );

    *x_index = -1;
    *y_index = -1;
    for_less( c, 0, N_DIMENSIONS )
    {
        if( x_axis[c] != 0.0 )
        {
            if( *x_index != -1 )
                return( FALSE );
            *x_index = c;
        }
        if( y_axis[c] != 0.0 )
        {
            if( *y_index != -1 )
                return( FALSE );
            *y_index = c;
        }
    }

    if( *x_index == *y_index )
        return( FALSE );

    *axis_index = N_DIMENSIONS - *x_index - *y_index;

    return( TRUE );
}
