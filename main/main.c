#include  <def_files.h>
#include  <def_string.h>
#include  <def_graphics.h>
#include  <def_globals.h>
#include  <def_arguments.h>

const  char  HARD_CODED_DISPLAY_DIRECTORY[] = "/usr/local/lib";
const  char  GLOBALS_FILENAME[]             = "display.globals";
const  char  MENU_FILENAME[]                = "display.menu";

int  main( argc, argv )
    int     argc;
    char    *argv[];
{
    char             *filename;
    String           runtime_directory;
    void             extract_directory();
    graphics_struct  *graphics;
    graphics_struct  *menu;
    Status           status;
    Status           initialize_graphics();
    Status           initialize_globals();
    Status           initialize_menu();
    Status           load_graphics_file();
    Status           create_graphics_window();
    Status           main_event_loop();
    Status           terminate_graphics();
    void             reset_view_parameters();
    void             update_view();
    void             set_model_scale();
    void             rebuild_selected_list();
    void             set_update_required();
    void             output_alloc_to_file();
    Status           delete_marching_cubes_table();
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

    status = initialize_globals( runtime_directory,
                                 HARD_CODED_DISPLAY_DIRECTORY,
                                 GLOBALS_FILENAME );

    if( status == OK )
    {
        status = initialize_graphics();
    }

    if( status == OK )
    {
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
        graphics->associated[SLICE_WINDOW] = (graphics_struct *) 0;

        menu->associated[THREE_D_WINDOW] = graphics;
        menu->associated[MENU_WINDOW] = menu;
        menu->associated[SLICE_WINDOW] = (graphics_struct *) 0;

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
        }
    }

    if( status == OK )
    {
        rebuild_selected_list( graphics, menu );
    }

    if( status == OK )
    {
        reset_view_parameters( graphics, &Default_line_of_sight,
                               &Default_horizontal );

        update_view( graphics );

        set_update_required( graphics, NORMAL_PLANES );
    }

    if( status == OK )
    {
        status = main_event_loop();
    }

    (void) terminate_graphics();

    if( status == OK )
    {
        status = delete_marching_cubes_table();
    }

    output_alloc_to_file( ".alloc_stats" );

    if( status != OK )
    {
        PRINT( "Program ended with error %d\n", (int) status );
    }

    return( (int) status );
}
