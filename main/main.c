#include  <def_display.h>

#define  HARD_CODED_DISPLAY_DIRECTORY     "/usr/local/lib"
#define  GLOBALS_FILENAME                 "display.globals"
#define  MENU_FILENAME                    "display.menu"


#define  GLOBALS_LOOKUP_NAME  display_globals
#include  <def_globals.h>
#define   REGISTER_GLOBALS_FILENAME   "display.globals"

public  Status  change_global_variable(
    char   str[],
    char   variable_name[],
    char   new_value[] )
{
    return( set_or_get_global_variable(
                   SIZEOF_STATIC_ARRAY(display_globals),
                   display_globals, str, variable_name, new_value ) );
}


int  main( argc, argv )
    int     argc;
    char    *argv[];
{
    char             *filename;
    String           runtime_directory;
    display_struct   *graphics;
    display_struct   *menu;
    Status           status;
    String           globals_filename;
    char             *title;

    if( argc == 1 )
        title = argv[0];
    else
        title = argv[1];

    if( getenv("DISPLAY_DIRECTORY") != (char *) 0 )
    {
        (void) strcpy( runtime_directory, getenv("DISPLAY_DIRECTORY") );
    }
    else
    {
        extract_directory( argv[0], runtime_directory );
    }

    get_absolute_filename( REGISTER_GLOBALS_FILENAME, runtime_directory,
                           globals_filename );

    status = OK;

    if( file_exists( globals_filename ) )
    {
        status = input_globals_file( SIZEOF_STATIC_ARRAY(display_globals),
                                     display_globals, globals_filename );
    }

    if( status == OK )
    {
        initialize_graphics();

        status = create_graphics_window( THREE_D_WINDOW,
                                         &graphics, title, 0, 0 );
    }

    if( status == OK )
    {
        status = create_graphics_window( MENU_WINDOW, &menu, title,
                                         Menu_window_width,
                                         Menu_window_height );
    }

    if( status == OK )
    {
        graphics->associated[THREE_D_WINDOW] = graphics;
        graphics->associated[MENU_WINDOW] = menu;
        graphics->associated[SLICE_WINDOW] = (display_struct *) 0;

        menu->associated[THREE_D_WINDOW] = graphics;
        menu->associated[MENU_WINDOW] = menu;
        menu->associated[SLICE_WINDOW] = (display_struct *) 0;

        status = initialize_menu( menu, runtime_directory,
                                  HARD_CODED_DISPLAY_DIRECTORY,
                                  MENU_FILENAME );
    }

    if( status == OK )
    {
        initialize_argument_processing( argc, argv );

        while( get_string_argument( "", &filename ) )
        {
            status = load_graphics_file( graphics, filename );
            if( status != OK )
                print( "Error loading %s\n", filename );
        }
    }

    if( status == OK )
        rebuild_selected_list( graphics, menu );

    if( status == OK )
    {
        reset_view_parameters( graphics, &Default_line_of_sight,
                               &Default_horizontal );

        update_view( graphics );

        set_update_required( graphics, NORMAL_PLANES );
    }

    if( status == OK )
        status = main_event_loop();

    (void) terminate_graphics();

    if( status == OK )
        delete_marching_cubes_table();

    output_alloc_to_file( ".alloc_stats" );

    if( status != OK )
        print( "Program ended with error %d\n", (int) status );

    return( (int) status );
}
