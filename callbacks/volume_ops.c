
#include  <display.h>

private  void  change_current_slice_by_one(
    display_struct   *display,
    int              delta )
{
    display_struct   *slice_window;
    Volume           volume;
    Real             voxel[N_DIMENSIONS], separations[N_DIMENSIONS];
    int              sizes[N_DIMENSIONS], axis_index, volume_index;

    if( get_slice_window( display, &slice_window ) &&
        get_axis_index_under_mouse( display, &volume_index, &axis_index ) )
    {
        volume = get_nth_volume( slice_window, volume_index );
        get_volume_sizes( volume, sizes );
        get_volume_separations( volume, separations );

        if( separations[axis_index] < 0.0 )
            delta = -delta;

        get_current_voxel( slice_window, volume_index, voxel );

        voxel[axis_index] = (Real) ROUND( voxel[axis_index] + (Real) delta );

        if( voxel[axis_index] < 0.0 )
            voxel[axis_index] = 0.0;
        else if( voxel[axis_index] > (Real) sizes[axis_index]-1.0 )
            voxel[axis_index] = (Real) sizes[axis_index] - 1.0;

        if( set_current_voxel( slice_window, volume_index, voxel ))
        {
            if( update_cursor_from_voxel( slice_window ) )
            {
                set_update_required( get_three_d_window(slice_window),
                                     get_cursor_bitplanes() );
            }
        }
    }
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(move_slice_plus)
{
    change_current_slice_by_one( display, 1 );

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(move_slice_plus )
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(move_slice_minus)
{
    change_current_slice_by_one( display, -1 );

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(move_slice_minus )
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(toggle_slice_visibility)
{
    int              view_index, volume_index;
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 0 &&
        get_slice_view_index_under_mouse( slice_window, &view_index ) )
    {
        volume_index = get_current_volume_index( slice_window );
        set_slice_visibility( slice_window, volume_index, view_index,
                !get_slice_visibility(slice_window,volume_index,view_index) );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(toggle_slice_visibility )
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(toggle_cross_section_visibility)
{
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) )
    {
        set_volume_cross_section_visibility( display,
                                !get_volume_cross_section_visibility(display) );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(toggle_cross_section_visibility )
{
    return( slice_window_exists(display) );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(reset_current_slice_view)
{
    int              view_index;
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 0 &&
        get_slice_view_index_under_mouse( slice_window, &view_index ) )
    {
        reset_slice_view( slice_window, view_index );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(reset_current_slice_view )
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(colour_code_objects )
{
    object_struct           *object, *current_object;
    Volume                  volume;
    object_traverse_struct  object_traverse;

    if( get_current_object(display,&current_object) &&
        get_slice_window_volume( display, &volume ) )
    {
        initialize_object_traverse( &object_traverse, 1, &current_object );

        while( get_next_object_traverse(&object_traverse,&object) )
        {
            colour_code_an_object( display, object );
        }

        set_update_required( display, NORMAL_PLANES );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(colour_code_objects )
{
    return( get_n_volumes(display) > 0 &&
            current_object_exists(display) );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(create_3d_slice)
{
    display_struct   *slice_window;
    int              x_index, y_index, axis_index, view_index;
    Real             current_voxel[N_DIMENSIONS];
    object_struct    *object;

    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 0 &&
        get_slice_view_index_under_mouse( slice_window, &view_index ) &&
        slice_has_ortho_axes( slice_window,
                         get_current_volume_index(slice_window), view_index,
                              &x_index, &y_index, &axis_index ) )

    {
        get_current_voxel( slice_window,
                        get_current_volume_index(slice_window), current_voxel );
        object = create_3d_slice_quadmesh( get_volume(display), axis_index,
                                           current_voxel[axis_index] );

        colour_code_an_object( display, object );

        add_object_to_current_model( display, object );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(create_3d_slice)
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(resample_slice_window_volume)
{
    int              sizes[N_DIMENSIONS];
    int              new_nx, new_ny, new_nz;
    char             label[1000];
    display_struct   *slice_window;
    Volume           volume, resampled_volume;

    if( get_slice_window_volume( display, &volume ) &&
        get_slice_window( display, &slice_window ) )
    {
        get_volume_sizes( volume, sizes );

        print( "The original volume is %d by %d by %d.\n",
               sizes[X], sizes[Y], sizes[Z] );

        print( "Enter desired resampled size: " );

        if( input_int( stdin, &new_nx ) == OK &&
            input_int( stdin, &new_ny ) == OK &&
            input_int( stdin, &new_nz ) == OK &&
            (new_nx > 0 || new_ny > 0 || new_nz > 0) )
        {
            resampled_volume = smooth_resample_volume(
                                        volume, new_nx, new_ny, new_nz );

            (void) sprintf( label, "Subsampled to %d,%d,%d: %s",
                            new_nx, new_ny, new_nz,
                            get_volume_filename(slice_window,
                                   get_current_volume_index(slice_window) ) );

            add_slice_window_volume( slice_window, label, resampled_volume );
        }

        (void) input_newline( stdin );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(resample_slice_window_volume)
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(box_filter_slice_window_volume)
{
    char             ch;
    Real             x_width, y_width, z_width;
    Real             separations[N_DIMENSIONS];
    char             label[1000];
    display_struct   *slice_window;
    Volume           volume, resampled_volume;

    if( get_slice_window_volume( display, &volume ) &&
        get_slice_window( display, &slice_window ) )
    {
        get_volume_separations( volume, separations );

        print( "Enter box filter  x_width, y_width, z_width, v/w: " );

        if( input_real( stdin, &x_width ) == OK &&
            input_real( stdin, &y_width ) == OK &&
            input_real( stdin, &z_width ) == OK &&
            input_nonwhite_character( stdin, &ch ) == OK &&
            (x_width > 1.0 || y_width > 1.0 || z_width > 1.0) )
        {
            (void) sprintf( label, "Box Filtered at %g,%g,%g,%c: %s",
                            x_width, y_width, z_width, ch,
                            get_volume_filename(slice_window,
                                   get_current_volume_index(slice_window) ) );

            if( ch == 'w' )
            {
                x_width /= ABS( separations[X] );
                y_width /= ABS( separations[Y] );
                z_width /= ABS( separations[Z] );
            }

            resampled_volume = create_box_filtered_volume( volume,
                                        NC_BYTE, FALSE, 0.0, 0.0,
                                        x_width, y_width, z_width );

            add_slice_window_volume( slice_window, label, resampled_volume );
        }

        (void) input_newline( stdin );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(box_filter_slice_window_volume)
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(pick_slice_angle_point)
{
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 0 )
    {
        start_picking_slice_angle( slice_window );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(pick_slice_angle_point)
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION( rotate_slice_axes )
{
    initialize_rotating_slice( display );

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(rotate_slice_axes )
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(reset_slice_crop)
{
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) )
    {
        reset_crop_box_position( slice_window );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(reset_slice_crop)
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(toggle_slice_crop_visibility)
{
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) )
    {
        toggle_slice_crop_box_visibility( slice_window );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(toggle_slice_crop_visibility)
{
    display_struct   *slice_window;
    BOOLEAN          visible;

    if( get_slice_window( display, &slice_window ) )
        visible = slice_window->slice.crop.crop_visible;
    else
        visible = OFF;

    set_menu_text_on_off( menu_window, menu_entry, visible );

    return( slice_window_exists(display) );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(pick_crop_box_edge)
{
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 0 )
    {
        start_picking_crop_box( slice_window );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(pick_crop_box_edge)
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(set_crop_box_filename)
{
    display_struct   *slice_window;
    STRING           filename;

    if( get_slice_window( display, &slice_window ) )
    {
        print( "Enter crop filename: " );

        if( input_string( stdin, filename, MAX_STRING_LENGTH, ' ' ) == OK )
        {
            set_crop_filename( slice_window, filename );
        }

        (void) input_newline( stdin );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(set_crop_box_filename)
{
    return( TRUE );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(load_cropped_volume)
{
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) )
    {
        crop_and_load_volume( slice_window );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(load_cropped_volume)
{
    return( slice_window_exists(display) );
}

private  void  do_histogram(
    display_struct   *display,
    BOOLEAN          labeled )
{
    int              x_index, y_index, view_index, axis_index;
    Real             voxel[MAX_DIMENSIONS], slice;
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 0 )
    {
        if( get_slice_view_index_under_mouse( slice_window, &view_index ) &&
            slice_has_ortho_axes( slice_window,
                       get_current_volume_index( slice_window ), view_index,
                       &x_index, &y_index, &axis_index ) )
        {
            get_current_voxel( slice_window,
                               get_current_volume_index(slice_window), voxel );
            slice = voxel[axis_index];
        }
        else
        {
            axis_index = -1;
            slice = 0.0;
        }

        compute_histogram( slice_window, axis_index, ROUND(slice), labeled );

        set_slice_viewport_update( slice_window, COLOUR_BAR_MODEL );
    }
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(redo_histogram)
{
    do_histogram( display, FALSE );

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(redo_histogram)
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(redo_histogram_labeled)
{
    do_histogram( display, TRUE );

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(redo_histogram_labeled)
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(print_voxel_origin)
{
    Real             voxel[MAX_DIMENSIONS], xw, yw, zw;
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 0 )
    {
        get_current_voxel( slice_window,
                      get_current_volume_index(slice_window), voxel );

        convert_voxel_to_world( get_volume(slice_window), voxel,
                                &xw, &yw, &zw );

        print( "Current world origin: %g %g %g\n", xw, yw, zw );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(print_voxel_origin)
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(print_slice_plane)
{
    display_struct   *slice_window;
    Vector           normal;
    int              view_index;
    Real             perp_axis[MAX_DIMENSIONS], xw, yw, zw;

    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 0 &&
        get_slice_view_index_under_mouse( slice_window, &view_index ) )
    {
        get_slice_perp_axis( slice_window,
                get_current_volume_index(slice_window), view_index, perp_axis );
        convert_voxel_vector_to_world( get_volume(slice_window),
                                       perp_axis, &xw, &yw, &zw );

        fill_Vector( normal, xw, yw, zw );
        NORMALIZE_VECTOR( normal, normal );

        print( "View: %d    Perpendicular: %g %g %g\n",
               view_index,
               Vector_x(normal), Vector_y(normal), Vector_z(normal) );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(print_slice_plane)
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(type_in_voxel_origin)
{
    Real             voxel[MAX_DIMENSIONS], xw, yw, zw;
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 0 )
    {
        print( "Enter x y z world coordinate: " );

        if( input_real( stdin, &xw ) == OK &&
            input_real( stdin, &yw ) == OK &&
            input_real( stdin, &zw ) == OK )
        {
            convert_world_to_voxel( get_volume(slice_window), xw, yw, zw,
                                    voxel );

            if( set_current_voxel( slice_window,
                        get_current_volume_index(slice_window), voxel ) )
            {
                if( update_cursor_from_voxel( slice_window ) )
                    set_update_required( display, NORMAL_PLANES );
            }
        }

        (void) input_newline( stdin );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(type_in_voxel_origin)
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(type_in_slice_plane)
{
    int              view_index;
    Real             perp_axis[MAX_DIMENSIONS], xw, yw, zw;
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 0 &&
        get_slice_view_index_under_mouse( slice_window, &view_index ) )
    {
        print( "View %d:  enter x y z plane normal in world coordinate: ",
               view_index );

        if( input_real( stdin, &xw ) == OK &&
            input_real( stdin, &yw ) == OK &&
            input_real( stdin, &zw ) == OK )
        {
            convert_world_vector_to_voxel( get_volume(slice_window), xw, yw, zw,
                                    perp_axis );

            set_slice_plane_perp_axis( slice_window,
                                       get_current_volume_index(slice_window),
                                       view_index, perp_axis);
            reset_slice_view( slice_window, view_index );
        }

        (void) input_newline( stdin );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(type_in_slice_plane)
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(toggle_slice_cross_section_visibility)
{
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) )
    {
        slice_window->slice.cross_section_visibility =
                             !slice_window->slice.cross_section_visibility;
        set_slice_cross_section_update( slice_window, -1 );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(toggle_slice_cross_section_visibility)
{
    return( slice_window_exists(display) );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(set_current_arbitrary_view)
{
    int              view_index;
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 0 &&
        get_slice_view_index_under_mouse( slice_window, &view_index ) )
    {
        print( "Current arbitrary view is now: %d\n", view_index );

        slice_window->slice.cross_section_index = view_index;
        set_slice_cross_section_update( slice_window, -1 );
        rebuild_volume_cross_section( slice_window );
        set_update_required( display, NORMAL_PLANES );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(set_current_arbitrary_view)
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(toggle_slice_anchor)
{
    int              c, view_index;
    Vector           axis;
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 0 )
    {
        if( slice_window->slice.cross_section_vector_present ||
            !get_slice_view_index_under_mouse( slice_window, &view_index ) ||
            view_index == get_arbitrary_view_index( slice_window ) )
        {
            slice_window->slice.cross_section_vector_present = FALSE;
            print( "Slice anchor turned off\n" );
        }
        else
        {
            slice_window->slice.cross_section_vector_present = TRUE;
            print( "Setting slice anchor to view index: %d\n", view_index );

            get_slice_cross_section_direction( slice_window, view_index,
                           get_arbitrary_view_index( slice_window ),
                           &axis );

            for_less( c, 0, N_DIMENSIONS )
                slice_window->slice.cross_section_vector[c] =
                                                     Vector_coord(axis,c);
        }
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(toggle_slice_anchor)
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(delete_current_volume)
{
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 0 )
    {
        delete_slice_window_volume( slice_window,
                                    get_current_volume_index(slice_window) );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(delete_current_volume)
{
    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(toggle_current_volume)
{
    int              current;
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 1 )
    {
        current = get_current_volume_index( slice_window );
        current = (current + 1) % slice_window->slice.n_volumes;
        set_current_volume_index( slice_window, current );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(toggle_current_volume)
{
    int              current_index;

    if( get_n_volumes(display) > 1 )
        current_index = get_current_volume_index( display ) + 1;
    else
        current_index = 1;

    set_menu_text_int( menu_window, menu_entry, current_index );

    return( get_n_volumes(display) > 1 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(set_current_volume_opacity)
{
    int              current;
    Real             opacity;
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 0 )
    {
        current = get_current_volume_index( slice_window );

        if( current >= 0 )
        {
            print( "Enter volume opacity ( 0.0 <= o <= 1.0 ): " );

            if( input_real( stdin, &opacity ) == OK &&
                opacity >= 0.0 )
            {
                set_volume_opacity( slice_window, current, opacity );
            }

            (void) input_newline( stdin );
        }
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(set_current_volume_opacity)
{
    Real             opacity;
    int              current_index;
    display_struct   *slice_window;

    current_index = get_current_volume_index( display );

    if( current_index >= 0 && get_slice_window( display, &slice_window ) )
        opacity = slice_window->slice.volumes[current_index].opacity;
    else
        opacity = 1.0;

    set_menu_text_real( menu_window, menu_entry, opacity );

    return( get_n_volumes(display) > 1 );
}

private  int  get_current_visible_volume(
    display_struct   *display )
{
    int              current_visible, volume_index, view;
    display_struct   *slice_window;

    current_visible = -1;
    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 0 )
    {
        for_less( view, 0, N_SLICE_VIEWS )
        {
            for_less( volume_index, 0, get_n_volumes(slice_window) )
            {
                if( get_slice_visibility( slice_window, volume_index, view ) )
                {
                    current_visible = volume_index;
                    break;
                }
            }
            if( current_visible >= 0 )
                break;
        }

        if( volume_index < get_n_volumes(slice_window) )
            current_visible = volume_index;
    }

    return( current_visible );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(next_volume_visible)
{
    int              current, view, volume_index;
    display_struct   *slice_window;
    BOOLEAN          all_invisible;

    if( get_slice_window( display, &slice_window ) &&
        get_n_volumes(slice_window) > 0 )
    {
        current = get_current_visible_volume( slice_window );
        current = (current + 1) % get_n_volumes(slice_window);

        for_less( view, 0, N_SLICE_VIEWS )
        {
            all_invisible = TRUE;
            for_less( volume_index, 0, get_n_volumes(slice_window) )
            {
                if( get_slice_visibility( slice_window, volume_index, view ) )
                    all_invisible = FALSE;
                set_slice_visibility( slice_window, volume_index, view, OFF );
            }

            if( !all_invisible || get_n_volumes(slice_window) == 1 )
                set_slice_visibility( slice_window, current, view, ON );
        }

        set_current_volume_index( slice_window, current );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(next_volume_visible)
{
    set_menu_text_int( menu_window, menu_entry,
                       get_current_visible_volume(display) + 1 );

    return( get_n_volumes(display) > 0 );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION(toggle_slice_interpolation)
{
    int              view_index, volume_index, continuity;
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) )
    {
        continuity = slice_window->slice.degrees_continuity;
        ++continuity;
        if( continuity == 1 )
            continuity = 2;
        else if( continuity == 3 )
            continuity = -1;

        slice_window->slice.degrees_continuity = continuity;

        set_slice_window_all_update( slice_window, -1, UPDATE_SLICE );
    }

    return( OK );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(toggle_slice_interpolation )
{
    int              continuity;
    display_struct   *slice_window;
    char             *name;

    if( get_slice_window( display, &slice_window ) )
        continuity = slice_window->slice.degrees_continuity;
    else
        continuity = Initial_slice_continuity;

    switch( continuity )
    {
    case 0:   name = "trilinear";   break;
    case 2:   name = "tricubic";   break;

    case -1:
    default:  name = "near neigh";   break;
    }

    set_menu_text( menu_window, menu_entry, name );

    return( TRUE );
}

