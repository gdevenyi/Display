
#include  <display.h>

#define    TOO_MANY_COLOURS_FOR_LABELS  10000

private  void   set_colour_of_label(
    display_struct    *slice_window,
    int               label,
    Colour            colour );

public  void  initialize_slice_colour_coding(
    display_struct    *slice_window )
{
    int   label;

    initialize_colour_coding( &slice_window->slice.colour_coding,
                              (Colour_coding_types) Initial_colour_coding_type,
                              Colour_below, Colour_above,
                              0.0, 1.0 );

    initialize_colour_bar( slice_window );

    slice_window->slice.label_colour_ratio = Label_colour_display_ratio;

    for_less( label, 0, NUM_LABELS )
    {
        slice_window->slice.colour_tables[label] = (Colour *) 0;
        slice_window->slice.label_colours_set[label] = FALSE;
    }
}

public  void  set_colour_coding_for_new_volume(
    display_struct    *slice_window )
{
    Real             low_limit, high_limit;
    Real             min_voxel, max_voxel;
    Real             min_value, max_value;
    int              label, ind;
    Volume           volume;
    BOOLEAN          too_many;
    Colour           col;
    static Colour    default_colours[] = { RED, GREEN, BLUE,
                                           CYAN, MAGENTA, YELLOW,
                                           BLUE_VIOLET, DEEP_PINK,
                                           GREEN_YELLOW, LIGHT_SEA_GREEN,
                                           MEDIUM_TURQUOISE, PURPLE };

    (void) get_slice_window_volume( slice_window, &volume );

    for_less( label, 0, NUM_LABELS )
    {
        slice_window->slice.colour_tables[label] = (Colour *) 0;
        slice_window->slice.label_colours_set[label] = FALSE;
    }

    get_volume_voxel_range( volume, &min_voxel, &max_voxel );

    too_many = ((int) max_voxel - (int) min_voxel + 1 >
                TOO_MANY_COLOURS_FOR_LABELS);

    if( too_many )
        slice_window->slice.n_labels = 1;
    else
        slice_window->slice.n_labels = NUM_LABELS;

    set_colour_of_label( slice_window, 0, WHITE );

    if( too_many )
    {
        for_less( label, 1, NUM_LABELS )
        {
            slice_window->slice.colour_tables[label] = 
                 slice_window->slice.colour_tables[0];
        }
    }
    else
    {
        for_less( label, 1, NUM_LABELS )
        {
            ind = (label - 1) % SIZEOF_STATIC_ARRAY(default_colours);

            col = default_colours[ind];
            if( label & get_active_bit() )
                col = SCALE_COLOUR( col, 0.5 );

            set_colour_of_label( slice_window, label, col );
        }

        set_colour_of_label( slice_window, get_label_bit(),
                             Labeled_voxel_colour );
    }

    get_volume_real_range( volume, &min_value, &max_value );

    rebuild_colour_tables( slice_window );

    low_limit = min_value + Initial_low_limit_position *
                (max_value - min_value);
    high_limit = min_value + Initial_high_limit_position *
                 (max_value - min_value);
    change_colour_coding_range( slice_window, low_limit, high_limit );

    rebuild_colour_bar( slice_window );
}

public  void  delete_slice_colour_coding(
    slice_window_struct   *slice )
{
    int      i;
    Real     min_voxel, max_voxel;
    Colour   *ptr;

    get_volume_voxel_range( slice->volume, &min_voxel, &max_voxel );

    for_less( i, 0, slice->n_labels )
    {
        if( slice->colour_tables[i] != (Colour *) 0 )
        {
            ptr = slice->colour_tables[i];
            ptr += (int) min_voxel;
            FREE( ptr );
        }
    }
}

public  void  change_colour_coding_range(
    display_struct    *slice_window,
    Real              min_value,
    Real              max_value )
{
    set_colour_coding_min_max( &slice_window->slice.colour_coding,
                               min_value, max_value );

    colour_coding_has_changed( slice_window );
}

public  void  create_colour_tables(
    display_struct    *slice_window,
    int               label )
{
    Real        min_voxel, max_voxel;
    Colour      *ptr;

    get_volume_voxel_range( slice_window->slice.volume,
                            &min_voxel, &max_voxel );

    ALLOC( ptr, (int) max_voxel - (int) min_voxel + 1 );

    slice_window->slice.colour_tables[label] = ptr - (int) min_voxel;
}

private  void   set_colour_of_label(
    display_struct    *slice_window,
    int               label,
    Colour            colour )
{
    slice_window->slice.label_colours[label] = colour;

    if( slice_window->slice.colour_tables[label] == (Colour *) 0 )
        create_colour_tables( slice_window, label );

    rebuild_colour_table_for_label( slice_window, label );
}

public  void   add_label_colour(
    display_struct    *slice_window,
    int               label,
    Colour            colour )
{
    set_colour_of_label( slice_window, label, colour );
    slice_window->slice.label_colours_set[label] = TRUE;
}

public  int  lookup_label_colour(
    display_struct    *slice_window,
    Colour            colour )
{
    BOOLEAN   found_colour, found_empty;
    int       i, first_empty, label;

    found_colour = FALSE;
    found_empty = FALSE;

    for_less( i, 1, MIN( get_active_bit(), slice_window->slice.n_labels ) )
    {
        label = i;

        if( slice_window->slice.label_colours != (Colour *) NULL )
        {
            if( equal_colours( slice_window->slice.label_colours[label],
                               colour ) )
            {
                found_colour = TRUE;
                break;
            }
        }

 
        if( !found_empty && !slice_window->slice.label_colours_set[label] )
        {
            found_empty = TRUE;
            first_empty = label;
        }
    }

    if( !found_colour )
    {
        if( found_empty )
        {
            label = first_empty;
            add_new_label( slice_window, label, colour );
        }
        else
            label = get_label_bit();
    }

    return( label );
}

private  Colour  apply_label_colour(
    display_struct    *slice_window,
    Colour            col,
    int               label )
{
    Colour           label_col, mult, scaled_col;

    if( label != 0 )
    {
        label_col = slice_window->slice.label_colours[label];
        MULT_COLOURS( mult, label_col, col );
        mult = SCALE_COLOUR( mult, 1.0-slice_window->slice.label_colour_ratio);
        scaled_col = SCALE_COLOUR( label_col,
                                   slice_window->slice.label_colour_ratio);
        ADD_COLOURS( col, mult, scaled_col );
    }

    return( col );
}

private  Colour  get_slice_colour_coding(
    display_struct    *slice_window,
    Real              value,
    int               label )
{
    Colour           col;

    col = get_colour_code( &slice_window->slice.colour_coding, value );

    if( label > 0 )
        col = apply_label_colour( slice_window, col, label );

    return( col );
}

public  void  rebuild_colour_table_for_label(
    display_struct    *slice_window,
    int               label )
{
    int              voxel;
    Real             value;
    Colour           colour;
    Real             min_voxel, max_voxel;

    if( label >= slice_window->slice.n_labels )
        return;

    get_volume_voxel_range( get_volume(slice_window), &min_voxel, &max_voxel );

    for_inclusive( voxel, (int) min_voxel, (int) max_voxel )
    {
        value = CONVERT_VOXEL_TO_VALUE( get_volume(slice_window), voxel );
        colour = get_slice_colour_coding( slice_window, value, label );

        slice_window->slice.colour_tables[label][voxel] = colour;
    }
}

public  void  rebuild_colour_tables(
    display_struct    *slice_window )
{
    int              label;

    for_less( label, 0, slice_window->slice.n_labels )
    {
        if( slice_window->slice.colour_tables[label] != (Colour *) 0 )
            rebuild_colour_table_for_label( slice_window, label );
    }
}

public  void  colour_coding_has_changed(
    display_struct    *display )
{
    display_struct    *slice_window;

    slice_window = display->associated[SLICE_WINDOW];

    if( slice_window != (display_struct  *) 0 )
    {
        rebuild_colour_tables( slice_window );

        rebuild_colour_bar( slice_window );
        set_slice_window_all_update( slice_window );
    }
}
