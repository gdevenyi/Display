#ifndef  DEF_GRAPHICS
#define  DEF_GRAPHICS

#include  <def_graphics_dependent.h>

#define  OFF  0
#define  ON   1

#define  PI  3.141592653589793

#define  DEG_TO_RAD   (PI / 180.0)
#define  RAD_TO_DEG   (180.0 / PI)

typedef  struct
{
    Boolean     perspective_flag;
    Point       origin;
    Vector      line_of_sight;
    Vector      horizontal;
    Vector      up;
    Real        front_distance, back_distance;
    Real        perspective_distance;
    Real        window_width, window_height;
} view_struct;

typedef enum { AMBIENT_LIGHT,
               DIRECTIONAL_LIGHT,
               POINT_LIGHT,
               SPOT_LIGHT
             } light_types;

typedef  struct
{
    light_types   light_type;
    Colour        colour;
    Vector        direction;
    Vector        x_axis;
    Vector        y_axis;
    Point         position;
    Real          spot_exponent;
    Real          spot_angle;
} light_struct;

typedef  enum  { WIREFRAME_MODE, SHADED_MODE } render_modes;

typedef  enum  { FLAT_SHADING, GOURAUD_SHADING, PHONG_SHADING } shading_types;

typedef  struct
{
    render_modes    render_mode;
    shading_types   shading_type;
    Boolean         master_light_switch;
}  render_struct;

typedef  struct
{
    Colour    colour;
    Surfprop  surfprop;

    int       n_points;
    Point     *points;
    Vector    *normals;

    int       n_triangles;
    int       *triangle_indices;

} triangles_struct;

#endif
