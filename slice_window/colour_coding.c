
#include  <display.h>

#define    MAX_LABEL_COLOUR_TABLE_SIZE    2000000

private  void  rebuild_colour_table(
    display_struct    *slice_window,
    int               volume_index );

private  BOOLEAN  is_shared_label_volume(
    slice_window_struct   *slice,
    int                   volume_index )
{
    int     i;

    for_less( i, 0, slice->n_volumes )
    {
        if( i != volume_index &&
            slice->volumes[i].labels == slice->volumes[volume_index].labels )
        {
            return( TRUE );
        }
    }

    return( FALSE );
}

private  void  delete_slice_labels(
    slice_window_struct   *slice,
    int                   volume_index )
{
    if( slice->volumes[volume_index].labels == (Volume) NULL )
        return;

    if( !is_shared_label_volume( slice, volume_index ) )
        delete_volume( slice->volumes[volume_index].labels );

    slice->volumes[volume_index].labels = NULL;
    slice->volumes[volume_index].labels_filename[0] = (char) 0;

    FREE( slice->volumes[volume_index].label_colour_table );
}

public  void  delete_slice_colour_coding(
    slice_window_struct   *slice,
    int                   volume_index )
{
    Colour   *ptr;

    if( slice->volumes[volume_index].labels == (Volume) NULL )
        return;

    delete_slice_labels( slice, volume_index );

    ptr = slice->volumes[volume_index].colour_table;

    if( ptr != NULL )
    {
        ptr += (int) slice->volumes[volume_index].offset;
        FREE( ptr );
    }
}

private  void  realloc_label_colour_table(
    display_struct    *slice_window,
    int               volume_index )
{
    int              ind;
    int              label, n_labels;
    static Colour    default_colours[] = { RED, GREEN, BLUE,
                                           CYAN, MAGENTA, YELLOW,
                                           BLUE_VIOLET, DEEP_PINK,
                                           GREEN_YELLOW, LIGHT_SEA_GREEN,
                                           MEDIUM_TURQUOISE, PURPLE };

    n_labels = get_num_labels( slice_window, volume_index );
    ALLOC( slice_window->slice.volumes[volume_index].label_colour_table,
           n_labels );

    for_less( label, 1, n_labels )
    {
        ind = (label - 1) % SIZEOF_STATIC_ARRAY(default_colours);
        set_colour_of_label( slice_window, volume_index, label,
                             default_colours[ind] );
    }

    set_colour_of_label( slice_window, volume_index, 0,
                         make_rgba_Colour(0,0,0,0) );
}

private  BOOLEAN  find_similar_labels(
    display_struct    *slice_window,
    int               volume_index,
    int               *orig_index )
{
    int                 i, k, m;
    int                 sizes[MAX_DIMENSIONS], this_size[MAX_DIMENSIONS];
    BOOLEAN             same;
    Volume              this_volume, volume;
    General_transform   *this_gen_transform, *gen_transform;
    Transform           *this_transform, *transform;

    this_volume = get_nth_volume( slice_window, volume_index );
    get_volume_sizes( this_volume, this_size );
    this_gen_transform = get_voxel_to_world_transform( this_volume );
    if( get_transform_type(this_gen_transform) != LINEAR )
        return( FALSE );
    this_transform = get_linear_transform_ptr( this_gen_transform );

    for_less( i, 0, get_n_volumes(slice_window) )
    {
        if( i == volume_index )
            continue;

        volume = get_nth_volume( slice_window, i );
        get_volume_sizes( volume, sizes );

        if( sizes[X] != this_size[X] ||
            sizes[Y] != this_size[Y] ||
            sizes[Z] != this_size[Z] )
            continue;

        gen_transform = get_voxel_to_world_transform( volume );
        if( get_transform_type(gen_transform) != LINEAR )
            continue;

        transform = get_linear_transform_ptr( gen_transform );

        same = TRUE;
        for_less( k, 0, 4 )
        for_less( m, 0, 4 )
        {
            if( Transform_elem(*this_transform,k,m) !=
                Transform_elem(*transform,k,m) )
                same = FALSE;
        }

        if( same )
            break;
    }

    if( i < get_n_volumes(slice_window) )
        *orig_index = i;

    return( i < get_n_volumes(slice_window) );
}
    

private  void  create_colour_coding(
    display_struct    *slice_window,
    int               volume_index )
{
    int                   orig_index;
    nc_type               type;
    slice_window_struct   *slice;
    Volume                label_volume;

    slice = &slice_window->slice;

    if( slice->volumes[volume_index].n_labels <= (1 << 8) )
        type = NC_BYTE;
    else if( slice->volumes[volume_index].n_labels <= (1 << 16) )
        type = NC_SHORT;
    else
        type = NC_LONG;

    /*--- this will break if you change the number of labels for a shared
          volume.  Later I will rewrite this code to handle this case */

    if( slice->share_labels_flag &&
        find_similar_labels( slice_window, volume_index, &orig_index ) )
    {
        slice->volumes[volume_index].labels = slice->volumes[orig_index].labels;
        (void) strcpy( slice->volumes[volume_index].labels_filename,
                       slice->volumes[orig_index].labels_filename );
    }
    else
    {
        slice->volumes[volume_index].labels = create_label_volume(
                          get_nth_volume( slice_window, volume_index ), type );
        slice->volumes[volume_index].labels_filename[0] = (char) 0;

        set_volume_voxel_range( slice->volumes[volume_index].labels, 0.0,
                            (Real) slice->volumes[volume_index].n_labels-1.0 );
    }

    realloc_label_colour_table( slice_window, volume_index );
}

public  void  set_slice_window_number_labels(
    display_struct    *slice_window,
    int               volume_index,
    int               n_labels )
{
    delete_slice_labels( &slice_window->slice, volume_index );

    slice_window->slice.volumes[volume_index].n_labels = n_labels;

    create_colour_coding( slice_window, volume_index );
}

private  void  alloc_colour_table(
    display_struct    *slice_window,
    int               volume_index )
{
    Real        min_voxel, max_voxel;
    Colour      *ptr;

    if( is_an_rgb_volume(get_nth_volume(slice_window,volume_index)) )
    {
        slice_window->slice.volumes[volume_index].colour_table = NULL;
        return;
    }

    get_volume_voxel_range( get_nth_volume(slice_window,volume_index),
                            &min_voxel, &max_voxel );

    ALLOC( ptr, (int) max_voxel - (int) min_voxel + 1 );

    slice_window->slice.volumes[volume_index].offset = (int) min_voxel;
    slice_window->slice.volumes[volume_index].colour_table =
                                              ptr - (int) min_voxel;
}

public  void  initialize_slice_colour_coding(
    display_struct    *slice_window,
    int               volume_index )
{
    Real             low_limit, high_limit;
    Real             min_value, max_value;

    initialize_colour_coding(
           &slice_window->slice.volumes[volume_index].colour_coding,
           (Colour_coding_types) Initial_colour_coding_type,
           Colour_below, Colour_above, 0.0, 1.0 );

    slice_window->slice.volumes[volume_index].label_colour_opacity =
                                                        Label_colour_opacity;
    slice_window->slice.volumes[volume_index].n_labels = Initial_num_labels;

    slice_window->slice.volumes[volume_index].offset = 0;
    slice_window->slice.volumes[volume_index].colour_table = (Colour *) NULL;
    slice_window->slice.volumes[volume_index].label_colour_table =
                                                 (Colour *) NULL;
    slice_window->slice.volumes[volume_index].labels = (Volume) NULL;
    slice_window->slice.volumes[volume_index].labels_filename[0] = (char) 0;

    alloc_colour_table( slice_window, volume_index );
    rebuild_colour_table( slice_window, volume_index );
    create_colour_coding( slice_window, volume_index );

    get_volume_real_range( get_nth_volume(slice_window,volume_index),
                           &min_value, &max_value );

    low_limit = min_value + Initial_low_limit_position *
                (max_value - min_value);
    high_limit = min_value + Initial_high_limit_position *
                 (max_value - min_value);

    change_colour_coding_range( slice_window, volume_index,
                                low_limit, high_limit );
}

public  Volume  get_nth_label_volume(
    display_struct   *display,
    int              volume_index )
{
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) )
    {
        return( slice_window->slice.volumes[volume_index].labels );
    }
    else
        return( (Volume) NULL );
}

public  Volume  get_label_volume(
    display_struct   *display )
{
    return( get_nth_label_volume( display,
                                  get_current_volume_index(display)) );
}

public  BOOLEAN  label_volume_exists(
    display_struct   *display )
{
    Volume   label;

    label = get_label_volume( display );

    return( label != NULL && label->data != NULL );
}

public  BOOLEAN  get_label_visibility(
    display_struct    *slice_window,
    int               volume_index,
    int               view_index )
{
    Volume   label;

    label = get_nth_label_volume( slice_window, volume_index );

    return( get_slice_visibility( slice_window, volume_index, view_index ) &&
            slice_window->slice.volumes[volume_index].display_labels &&
            label != NULL && label->data != NULL );
}

public  int  get_num_labels(
    display_struct   *display,
    int              volume_index )
{
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window ) == NULL ||
        slice_window->slice.n_volumes == 0 )
        return( Initial_num_labels );
    else
        return( slice_window->slice.volumes[volume_index].n_labels );
}

private  Colour  apply_label_colour(
    display_struct    *slice_window,
    int               volume_index,
    Colour            col,
    int               label )
{
    Real      r1, g1, b1, a1, r2, g2, b2, a2;
    Real      r, g, b, a;
    Colour    label_col;

    if( label != 0 )
    {
        label_col = slice_window->slice.volumes[volume_index].
                                                 label_colour_table[label];

        r1 = get_Colour_r_0_1(col);
        g1 = get_Colour_g_0_1(col);
        b1 = get_Colour_b_0_1(col);
        a1 = get_Colour_a_0_1(col);

        r2 = get_Colour_r_0_1(label_col);
        g2 = get_Colour_g_0_1(label_col);
        b2 = get_Colour_b_0_1(label_col);
        a2 = get_Colour_a_0_1(label_col);

        r = a1 * r1 * (1.0 - a2) + a2 * r2;
        g = a1 * g1 * (1.0 - a2) + a2 * g2;
        b = a1 * b1 * (1.0 - a2) + a2 * b2;
        a = a1 * (1.0 - a2) + a2 * a2;

        col = make_rgba_Colour_0_1( r, g, b, a );
    }

    return( col );
}

private  Colour  get_slice_colour_coding(
    display_struct    *slice_window,
    int               volume_index,
    Real              value,
    int               label )
{
    Colour           col;

    col = get_colour_code( &slice_window->slice.volumes[volume_index].
                           colour_coding, value );

    if( label > 0 )
        col = apply_label_colour( slice_window, volume_index, col, label );

    return( col );
}

private  void  rebuild_colour_table(
    display_struct    *slice_window,
    int               volume_index )
{
    Volume           volume;
    int              voxel;
    Real             value, r, g, b, a, opacity;
    Colour           colour;
    Real             min_voxel, max_voxel;

    volume = get_nth_volume(slice_window,volume_index);

    if( is_an_rgb_volume(volume) )
        return;

    get_volume_voxel_range( volume, &min_voxel, &max_voxel );

    opacity = slice_window->slice.volumes[volume_index].opacity;

    for_inclusive( voxel, (int) min_voxel, (int) max_voxel )
    {
        value = CONVERT_VOXEL_TO_VALUE( volume, voxel );
        colour = get_colour_code( &slice_window->slice.volumes[volume_index].
                                  colour_coding, value );

        r = get_Colour_r_0_1( colour );
        g = get_Colour_g_0_1( colour );
        b = get_Colour_b_0_1( colour );
        a = get_Colour_a_0_1( colour );

        slice_window->slice.volumes[volume_index].colour_table[voxel] =
                    make_rgba_Colour_0_1( r, g, b, a * opacity );
    }
}

public  void   set_colour_of_label(
    display_struct    *slice_window,
    int               volume_index,
    int               label,
    Colour            colour )
{
    Real  r, g, b;

    if( get_Colour_a(colour) == 255 )
    {
        r = get_Colour_r_0_1( colour );
        g = get_Colour_g_0_1( colour );
        b = get_Colour_b_0_1( colour );

        colour = make_rgba_Colour_0_1( r, g, b,
              slice_window->slice.volumes[volume_index].label_colour_opacity );
    }

    slice_window->slice.volumes[volume_index].label_colour_table[label] =colour;
}

public  Colour   get_colour_of_label(
    display_struct    *slice_window,
    int               volume_index,
    int               label )
{
    return( slice_window->slice.volumes[volume_index].
                                             label_colour_table[label] );
}

public  void   set_volume_opacity(
    display_struct   *slice_window,
    int              volume_index,
    Real             opacity )
{
    slice_window->slice.volumes[volume_index].opacity = opacity;

    rebuild_colour_table( slice_window, volume_index );

    colour_coding_has_changed( slice_window, volume_index, UPDATE_BOTH );
}

public  void   set_label_opacity(
    display_struct   *slice_window,
    int              volume_index,
    Real             opacity )
{
    int     i, n_labels;
    Real    r, g, b;
    Colour  *table;

    slice_window->slice.volumes[volume_index].label_colour_opacity = opacity;

    n_labels = get_num_labels(slice_window,volume_index);
    table = slice_window->slice.volumes[volume_index].label_colour_table;

    for_less( i, 1, n_labels )
    {
        r = get_Colour_r_0_1( table[i] );
        g = get_Colour_g_0_1( table[i] );
        b = get_Colour_b_0_1( table[i] );
        table[i] = make_rgba_Colour_0_1( r, g, b, opacity );
    }

    colour_coding_has_changed( slice_window, volume_index, UPDATE_LABELS );
}

public  void  colour_coding_has_changed(
    display_struct    *display,
    int               volume_index,
    Update_types      type )
{
    display_struct    *slice_window;

    if( get_slice_window( display, &slice_window ) )
    {
        if( type == UPDATE_SLICE || type == UPDATE_BOTH )
        {
            rebuild_colour_table( slice_window, volume_index );

            if( volume_index == get_current_volume_index(slice_window) )
                set_colour_bar_update( slice_window );
        }

        set_slice_window_all_update( slice_window, volume_index, type );
    }
}

public  void  change_colour_coding_range(
    display_struct    *slice_window,
    int               volume_index,
    Real              min_value,
    Real              max_value )
{
    set_colour_coding_min_max( &slice_window->slice.volumes[volume_index].
                              colour_coding,
                               min_value, max_value );

    colour_coding_has_changed( slice_window, volume_index, UPDATE_SLICE );
}

private  void  colour_code_points(
    display_struct        *slice_window,
    int                   continuity,
    Colour_flags          *colour_flag,
    Colour                *colours[],
    int                   n_points,
    Point                 points[] )
{
    int      i, int_voxel[MAX_DIMENSIONS], label, volume_index;
    Real     val, voxel[MAX_DIMENSIONS];
    Volume   volume, label_volume;
    Colour   colour, volume_colour;

    if( *colour_flag != PER_VERTEX_COLOURS )
    {
        if( n_points > 0 )
        {
            REALLOC( *colours, n_points );
        }
        else
        {
            FREE( *colours );
        }
        *colour_flag = PER_VERTEX_COLOURS;
    }

    for_less( i, 0, n_points )
    {
        colour = make_rgba_Colour( 0, 0, 0, 0 );

        for_less( volume_index, 0, slice_window->slice.n_volumes )
        {
            if( slice_window->slice.volumes[volume_index].opacity == 0.0 )
                continue;

            volume = get_nth_volume( slice_window, volume_index );
            label_volume = get_nth_label_volume( slice_window, volume_index );

            if( is_an_rgb_volume( volume ) )
            {
                convert_real_to_int_voxel( get_volume_n_dimensions(volume),
                                           voxel, int_voxel );
                if( int_voxel_is_within_volume( volume, int_voxel ) )
                {
                    GET_VOXEL_3D( volume_colour,
                                  volume, int_voxel[0], int_voxel[1],
                                  int_voxel[2] );
                }
            }
            else
            {
                (void) evaluate_volume_in_world( volume,
                                          Point_x(points[i]),
                                          Point_y(points[i]),
                                          Point_z(points[i]), continuity,
                                          &val, (Real *) NULL,
                                          (Real *) NULL, (Real *) NULL,
                                          (Real *) NULL, (Real *) NULL,
                                          (Real *) NULL, (Real *) NULL,
                                          (Real *) NULL, (Real *) NULL );
            }

            if( slice_window->slice.volumes[volume_index].display_labels &&
                label_volume != NULL && label_volume->data != NULL )
            {
                convert_world_to_voxel( volume,
                                        Point_x(points[i]),
                                        Point_y(points[i]),
                                        Point_z(points[i]), voxel );
                convert_real_to_int_voxel( get_volume_n_dimensions(volume),
                                           voxel, int_voxel );

                label = get_volume_label_data( label_volume, int_voxel );
            }
            else
                label = 0;

            if( is_an_rgb_volume( volume ) )
            {
                volume_colour = apply_label_colour( slice_window,
                                     volume_index, volume_colour, label );
            }
            else
            {
                volume_colour = get_slice_colour_coding( slice_window,
                                         volume_index, val, label );
            }

            COMPOSITE_COLOURS( colour, volume_colour, colour )
        }

        (*colours)[i] = colour;
    }
}

private  void  colour_code_object_points(
    display_struct         *slice_window,
    int                    continuity,
    object_struct          *object )
{
    Colour                  *colours;
    Colour_flags            colour_flag;
    polygons_struct         *polygons;
    quadmesh_struct         *quadmesh;
    lines_struct            *lines;
    marker_struct           *marker;

    switch( object->object_type )
    {
    case POLYGONS:
        polygons = get_polygons_ptr( object );
        colour_code_points( slice_window, continuity,
                            &polygons->colour_flag, &polygons->colours,
                            polygons->n_points, polygons->points );
        break;

    case QUADMESH:
        quadmesh = get_quadmesh_ptr( object );
        colour_code_points( slice_window, continuity,
                            &quadmesh->colour_flag, &quadmesh->colours,
                            quadmesh->m * quadmesh->n,
                            quadmesh->points );
        break;

    case LINES:
        lines = get_lines_ptr( object );
        colour_code_points( slice_window, continuity,
                            &lines->colour_flag, &lines->colours,
                            lines->n_points, lines->points );

    case MARKER:
        marker = get_marker_ptr( object );
        colour_flag = PER_VERTEX_COLOURS;
        colours = &marker->colour;
        colour_code_points( slice_window, continuity,
                            &colour_flag, &colours, 1, &marker->position );
        break;
    }
}

public  void  colour_code_an_object(
    display_struct   *display,
    object_struct    *object )
{
    display_struct   *slice_window;

    if( get_slice_window( display, &slice_window) )
        colour_code_object_points( slice_window, Volume_continuity, object );
}
