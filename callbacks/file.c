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
static char rcsid[] = "$Header: /private-cvsroot/visualization/Display/callbacks/file.c,v 1.18 1995-09-13 13:25:15 david Exp $";
#endif


#include  <display.h>

/* ARGSUSED */

public  DEF_MENU_FUNCTION( load_file )
{
    Status   status;
    STRING   filename;

    print( "Enter filename: " );

    status = input_string( stdin, filename, MAX_STRING_LENGTH, ' ' );

    (void) input_newline( stdin );

    if( status == OK )
        status = load_graphics_file( display, filename, FALSE );

    if( status == OK )
        graphics_models_have_changed( display );

    return( status );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(load_file )
{
    return( TRUE );
}

/* ARGSUSED */

public  DEF_MENU_FUNCTION( save_file )
{
    int            n_objects;
    object_struct  **object_list;
    object_struct  *current_object;
    Status         status;
    STRING         filename;

    status = OK;

    if( get_current_object( display, &current_object ) )
    {
        print( "Enter filename: " );

        status = input_string( stdin, filename, MAX_STRING_LENGTH, ' ' );

        (void) input_newline( stdin );

        if( status == OK && !check_clobber_file_default_suffix( filename,
                                                                "obj" ) )
            status = ERROR;

        if( current_object->object_type == MODEL )
        {
            n_objects = get_model_ptr(current_object)->n_objects;
            object_list = get_model_ptr(current_object)->objects;
        }
        else
        {
            n_objects = 1;
            object_list = &current_object;
        }

        if( status == OK )
        {
            status = output_graphics_file( filename, (File_formats) Save_format,
                                           n_objects, object_list );
            print( "Done saving.\n" );
        }
    }

    return( status );
}

/* ARGSUSED */

public  DEF_MENU_UPDATE(save_file )
{
    return( TRUE );
}
