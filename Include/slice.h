
#ifndef  DEF_SLICE
#define  DEF_SLICE

#include  <mni.h>
#include  <atlas.h>

typedef  struct
{
    int           axis_map[N_DIMENSIONS];
    Real          x_trans, y_trans;
    Real          x_scaling, y_scaling;
    Real          lower_view_limits[2];
    Real          upper_view_limits[2];
    BOOLEAN       update_flag;
    Filter_types  filter_type;
    Real          filter_width;
} slice_view_struct;

typedef  struct
{
    int   voxel_indices[3];
    int   id;
} label_struct;

typedef  enum  { FOUR_NEIGHBOURS, EIGHT_NEIGHBOURS } Neighbour_types;

typedef  struct
{
    int               n_labels;
    label_struct      *labels;
    int               min_threshold;
    int               max_threshold;
    Neighbour_types   connectivity;
} segmenting_struct;

typedef struct
{
    Real             top_offset;
    Real             bottom_offset;
    Real             left_offset;
    Real             bar_width;
    Real             tick_width;
    int              desired_n_intervals;
} colour_bar_struct;

#define  NUM_LABELS   256

typedef  struct
{
    Volume                 original_volume;
    Volume                 original_labels;

    Volume                 volume;
    Volume                 labels;

    Colour                 *colour_tables[NUM_LABELS];
    Real                   label_colour_ratio;
    BOOLEAN                label_colours_used[NUM_LABELS];
    Colour                 label_colours[NUM_LABELS];
    colour_coding_struct   colour_coding;
    colour_bar_struct      colour_bar;
    BOOLEAN                display_labels;

    int                    x_split, y_split;

    Real                   slice_index[N_DIMENSIONS];
    BOOLEAN                slice_locked[N_DIMENSIONS];
    slice_view_struct      slice_views[3];
    int                    next_to_update;

    segmenting_struct      segmenting;
    atlas_struct           atlas;

    Real                   x_brush_radius, y_brush_radius, z_brush_radius;
    int                    current_paint_label;
    object_struct          *brush_outline;

} slice_window_struct;


#endif
