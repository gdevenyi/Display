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
static char rcsid[] = "$Header: /private-cvsroot/visualization/Display/atlas/atlas.c,v 1.21 1996-07-02 12:56:11 david Exp $";
#endif

#include  <display.h>

private  atlas_position_struct  *get_closest_atlas_slice(
    int           axis,
    Real          slice_position,
    atlas_struct  *atlas );
private  Status  input_pixel_map(
    STRING       default_directory,
    STRING       image_filename,
    pixels_struct  *pixels );
private  BOOLEAN  find_appropriate_atlas_image(
    pixels_struct           *atlas_images,
    atlas_position_struct   *atlas_page,
    Real                    x_n_pixels,
    Real                    y_n_pixels,
    unsigned char           *atlas_image[],
    int                     *atlas_x_size,
    int                     *atlas_y_size,
    Real                    *x_atlas_to_voxel,
    Real                    *y_atlas_to_voxel );

private  const  int   ATLAS_SIZE[N_DIMENSIONS] = { 256, 256, 80 };
private  const  Real  ATLAS_STEPS[N_DIMENSIONS] = { 0.67, 0.86, 1.5 };
private  const  Real  ATLAS_STARTS[N_DIMENSIONS] = {  -86.095, -126.51, -37.5 };

public  void  initialize_atlas(
    atlas_struct   *atlas )
{
    atlas->input = FALSE;
    atlas->enabled = FALSE;
    atlas->opacity = Initial_atlas_opacity;
    atlas->transparent_threshold = Initial_atlas_transparent_threshold;
    atlas->slice_tolerance[X] = Initial_atlas_tolerance_x;
    atlas->slice_tolerance[Y] = Initial_atlas_tolerance_y;
    atlas->slice_tolerance[Z] = Initial_atlas_tolerance_z;
    atlas->flipped[X] = FALSE;
    atlas->flipped[Y] = FALSE;
    atlas->flipped[Z] = FALSE;
    atlas->slice_lookup[X] = (atlas_position_struct **) 0;
    atlas->slice_lookup[Y] = (atlas_position_struct **) 0;
    atlas->slice_lookup[Z] = (atlas_position_struct **) 0;
    atlas->n_pixel_maps = 0;
    atlas->n_pages = 0;
}

public  void  delete_atlas(
    atlas_struct   *atlas )
{
    int    axis, p;

    for_less( axis, 0, N_DIMENSIONS )
    {
        if( atlas->slice_lookup[axis] != (atlas_position_struct **) 0 )
            FREE( atlas->slice_lookup[axis] );
    }

    if( atlas->n_pixel_maps > 0 )
    {
        for_less( p, 0, atlas->n_pixel_maps )
            delete_pixels( &atlas->pixel_maps[p] );

        FREE( atlas->pixel_maps );
    }

    if( atlas->n_pages > 0 )
    {
        for_less( p, 0, atlas->n_pages )
        {
            if( atlas->pages[p].n_resolutions > 0 )
                FREE( atlas->pages[p].pixel_map_indices );
        }

        FREE( atlas->pages );
    }
}

private  Status  input_atlas(
    atlas_struct   *atlas,
    STRING         filename )
{
    Status           status;
    FILE             *file;
    STRING           *image_filenames, image_filename;
    char             axis_letter;
    Real             mm_position;
    int              pixel_index, page_index, axis_index;
    STRING           atlas_directory;
    progress_struct  progress;

    atlas_directory = extract_directory( filename );

    status = open_file( filename, READ_FILE, ASCII_FORMAT, &file );

    if( status == OK )
    {
        atlas->n_pixel_maps = 0;
        atlas->n_pages = 0;
        image_filenames = (STRING *) 0;

        while( input_string( file, &image_filename, ' ' ) == OK )
        {
            status = ERROR;

            if( input_nonwhite_character( file, &axis_letter ) != OK )
                break;

            if( axis_letter >= 'x' && axis_letter <= 'z' )
                axis_index = axis_letter - 'x';
            else if( axis_letter >= 'X' && axis_letter <= 'Z' )
                axis_index = axis_letter - 'X';
            else
                break;

            if( input_real( file, &mm_position ) != OK )
                break;

            for_less( pixel_index, 0, atlas->n_pixel_maps )
            {
                if( equal_strings(image_filenames[pixel_index],image_filename) )
                    break;
            }

            if( pixel_index == atlas->n_pixel_maps )
            {
                SET_ARRAY_SIZE( image_filenames,
                                atlas->n_pixel_maps, atlas->n_pixel_maps+1,
                                DEFAULT_CHUNK_SIZE );
                image_filenames[atlas->n_pixel_maps] = image_filename;
                SET_ARRAY_SIZE( atlas->pixel_maps,
                                atlas->n_pixel_maps, atlas->n_pixel_maps+1,
                                    DEFAULT_CHUNK_SIZE );

                ++atlas->n_pixel_maps;
            }

            for_less( page_index, 0, atlas->n_pages )
            {
                if( atlas->pages[page_index].axis == axis_index &&
                    atlas->pages[page_index].axis_position == mm_position )
                    break;
            }

            if( page_index == atlas->n_pages )
            {
                SET_ARRAY_SIZE( atlas->pages, atlas->n_pages,
                                atlas->n_pages+1, DEFAULT_CHUNK_SIZE );
                atlas->pages[atlas->n_pages].axis = axis_index;
                atlas->pages[atlas->n_pages].axis_position = mm_position;
                atlas->pages[atlas->n_pages].n_resolutions = 0;
                ++atlas->n_pages;
            }

            ADD_ELEMENT_TO_ARRAY(atlas->pages[page_index].pixel_map_indices,
                                 atlas->pages[page_index].n_resolutions,
                                 pixel_index, DEFAULT_CHUNK_SIZE );

            status = OK;
        }
    }

    if( status == OK )
        status = close_file( file );

    if( status == OK )
    {
        initialize_progress_report( &progress, FALSE, atlas->n_pixel_maps,
                                    "Reading Atlas" );

        for_less( pixel_index, 0, atlas->n_pixel_maps )
        {
            status = input_pixel_map( atlas_directory,
                                      image_filenames[pixel_index],
                                      &atlas->pixel_maps[pixel_index] );
            if( status != OK )
                break;

            update_progress_report( &progress, pixel_index+1 );
        }

        terminate_progress_report( &progress );
    }

    if( status == OK && atlas->n_pixel_maps > 0 )
    {
        for_less( pixel_index, 0, atlas->n_pixel_maps )
            delete_string( image_filenames[pixel_index] );

        FREE( image_filenames );
    }

    if( status == OK )
        atlas->input = TRUE;

    if( status != OK )
        print( "Error inputting atlas.\n" );

    delete_string( atlas_directory );

    return( status );
}

private  Status  input_pixel_map(
    STRING         default_directory,
    STRING         image_filename,
    pixels_struct  *pixels )
{
    Status         status;
    STRING         absolute_filename;
    File_formats   format;
    Object_types   object_type;
    BOOLEAN        eof;
    FILE           *file;

    absolute_filename = get_absolute_filename( image_filename,
                                               default_directory );

    status = open_file( absolute_filename, READ_FILE, BINARY_FORMAT, &file );

    if( status == OK )
        status = input_object_type( file, &object_type, &format, &eof );

    if( status == OK && !eof && object_type == PIXELS )
        status = io_pixels( file, READ_FILE, format, pixels );
    else
        status = ERROR;

    if( status == OK && pixels->pixel_type != COLOUR_INDEX_8BIT_PIXEL )
        status = ERROR;

    if( status == OK )
        status = close_file( file );

    delete_string( absolute_filename );

    return( status );
}

public  void  regenerate_atlas_lookup(
    display_struct    *slice_window )
{
    Volume            volume;
    atlas_struct      *atlas;
    Real              voxel[MAX_DIMENSIONS], world[N_DIMENSIONS];
    int               sizes[MAX_DIMENSIONS], axis, i, used_i;

    (void) get_slice_window_volume( slice_window, &volume );
    get_volume_sizes( volume, sizes );
    atlas = &slice_window->slice.atlas;

    for_less( axis, 0, N_DIMENSIONS )
    {
        if( atlas->slice_lookup[axis] != (atlas_position_struct **) 0 )
            FREE( atlas->slice_lookup[axis] );

        ALLOC( atlas->slice_lookup[axis], sizes[axis] );

        for_less( i, 0, sizes[axis] )
        {
            voxel[X] = 0.0;
            voxel[Y] = 0.0;
            voxel[Z] = 0.0;
            voxel[axis] = (Real) i;
            convert_voxel_to_world( volume, voxel,
                                    &world[X], &world[Y], &world[Z] );
            if( axis == X && atlas->flipped[X] )
                used_i = sizes[axis]-1-i;
            else
                used_i = i;

            atlas->slice_lookup[axis][used_i] = get_closest_atlas_slice(
                                                axis, world[axis], atlas );
        }
    }
}

private  Real  get_distance_from_voxel(
    Real  slice_position,
    Real  mm_coordinate )
{
    Real   distance;

    distance = FABS( slice_position - mm_coordinate );

    return( distance );
}

private  atlas_position_struct  *get_closest_atlas_slice(
    int           axis,
    Real          slice_position,
    atlas_struct  *atlas )
{
    int                     i;
    Real                    min_dist, dist;
    atlas_position_struct   *closest_so_far;

    closest_so_far = (atlas_position_struct *) 0;
    min_dist = 0.0;

    for_less( i, 0, atlas->n_pages )
    {
        if( atlas->pages[i].axis == axis )
        {
            dist = get_distance_from_voxel( slice_position,
                                            atlas->pages[i].axis_position );

            if( dist <= atlas->slice_tolerance[axis] &&
                (closest_so_far == NULL ||
                 dist < min_dist) )
            {
                closest_so_far = &atlas->pages[i];
                min_dist = dist;
            }
        }
    }

    return( closest_so_far );
}

public  void  set_atlas_state(
    display_struct    *slice_window,
    BOOLEAN           state )
{
    Status   status;

    status = OK;

    if( state && !slice_window->slice.atlas.input )
    {
        status = input_atlas( &slice_window->slice.atlas, Atlas_filename );

        if( status == OK )
            regenerate_atlas_lookup( slice_window );
    }

    slice_window->slice.atlas.enabled = state;
}

public  BOOLEAN  render_atlas_slice_to_pixels(
    atlas_struct  *atlas,
    Colour        image[],
    int           image_x_size,
    int           image_y_size,
    Real          voxel_start_indices[N_DIMENSIONS],
    int           a1,
    int           a2,
    int           axis_index,
    Real          x_pixel_to_voxel,
    Real          y_pixel_to_voxel,
    int           x_volume_size,
    int           y_volume_size )
{
    int            x, y, atlas_x_size, atlas_y_size;
    Real           x_atlas_to_voxel, y_atlas_to_voxel;
    int            x_pixel, y_pixel;
    Real           x_pixel_start, y_pixel_start;
    int            r_atlas, g_atlas, b_atlas;
    int            transparent_threshold;
    int            *x_pixels;
    Colour         atlas_pixel, *lookup, *pixels;
    unsigned  char *atlas_image;
    Real           opacity;

    if( !atlas->enabled || atlas->opacity <= 0.0 ||
        atlas->slice_lookup[axis_index] == (atlas_position_struct **) 0 ||
        atlas->slice_lookup[axis_index]
        [ROUND(voxel_start_indices[axis_index])] == (atlas_position_struct *) 0 ||
        !find_appropriate_atlas_image( atlas->pixel_maps,
                      atlas->slice_lookup[axis_index]
                  [ROUND(voxel_start_indices[axis_index])],
                  (Real) x_volume_size / x_pixel_to_voxel / ATLAS_STEPS[a1],
                  (Real) y_volume_size / y_pixel_to_voxel / ATLAS_STEPS[a2],
                  &atlas_image, &atlas_x_size, &atlas_y_size,
                  &x_atlas_to_voxel, &y_atlas_to_voxel ) )
    {
        return( FALSE );
    }

    opacity = atlas->opacity;
    transparent_threshold = atlas->transparent_threshold;
    lookup = get_8bit_rgb_pixel_lookup();

    x_pixel_start = voxel_start_indices[a1] / x_pixel_to_voxel /
                    x_atlas_to_voxel / ATLAS_STEPS[a1];
    y_pixel_start = voxel_start_indices[a2] / y_pixel_to_voxel /
                    y_atlas_to_voxel / ATLAS_STEPS[a2];

    ALLOC( x_pixels, image_x_size );

    for_less( x, 0, image_x_size )
    {
        x_pixels[x] = ROUND( x_pixel_start + (Real) x / x_atlas_to_voxel /
                                             ATLAS_STEPS[a1] );
        if( axis_index != X && atlas->flipped[axis_index] )
            x_pixels[x] = atlas_x_size - 1 - x_pixels[x];
    }

    for_less( y, 0, image_y_size )
    {
        pixels = &image[IJ(y,0,image_x_size)];

        y_pixel = ROUND( y_pixel_start + (Real) y / y_atlas_to_voxel /
                                         ATLAS_STEPS[a2] );

        if( y_pixel >= 0 && y_pixel < atlas_y_size )
        {
            for_less( x, 0, image_x_size )
            {
                x_pixel = x_pixels[x];

                if( x_pixel >= 0 && x_pixel < atlas_x_size )
                {
                    atlas_pixel = lookup[ (unsigned long) atlas_image
                                           [IJ(y_pixel,x_pixel,atlas_x_size)]];
                    r_atlas = get_Colour_r(atlas_pixel);
                    g_atlas = get_Colour_g(atlas_pixel);
                    b_atlas = get_Colour_b(atlas_pixel);

                    if( r_atlas <= transparent_threshold ||
                        g_atlas <= transparent_threshold ||
                        b_atlas <= transparent_threshold )
                    {
                        *pixels = make_rgba_Colour( r_atlas, g_atlas, b_atlas,
                                                    (int) (255.0 * opacity) );
                    }
                    else
                        *pixels = make_rgba_Colour( 0, 0, 0, 0 );
                }

                ++pixels;
            }
        }
    }

    FREE( x_pixels );

    return( TRUE );
}

private  BOOLEAN  find_appropriate_atlas_image(
    pixels_struct           *atlas_images,
    atlas_position_struct   *atlas_page,
    Real                    x_n_pixels,
    Real                    y_n_pixels,
    unsigned char           *atlas_image[],
    int                     *atlas_x_size,
    int                     *atlas_y_size,
    Real                    *x_atlas_to_voxel,
    Real                    *y_atlas_to_voxel )
{
    int             i, image_index;
    pixels_struct   *pixels;

    image_index = -1;

    for_less( i, 0, atlas_page->n_resolutions )
    {
        pixels = &atlas_images[atlas_page->pixel_map_indices[i]];

        if( pixels->x_size <= (int) x_n_pixels &&
            pixels->y_size <= (int) y_n_pixels &&
            (image_index < 0 || pixels->x_size > *atlas_x_size ||
             pixels->y_size > *atlas_y_size) )
        {
            image_index = i;
            *atlas_x_size = pixels->x_size;
            *atlas_y_size = pixels->y_size;
        }
    }

    if( image_index >= 0 )
    {
        *atlas_image =
              atlas_images[atlas_page->pixel_map_indices[image_index]].
                           data.pixels_8bit_colour_index;
        *x_atlas_to_voxel = (Real) x_n_pixels / (Real) *atlas_x_size;
        *y_atlas_to_voxel = (Real) y_n_pixels / (Real) *atlas_y_size;
    }

    return( image_index >= 0 );
}

public  BOOLEAN  is_atlas_loaded(
    display_struct  *display )
{
    display_struct  *slice_window;

    return( get_slice_window( display, &slice_window ) &&
            slice_window->slice.atlas.n_pixel_maps > 0 );
}
