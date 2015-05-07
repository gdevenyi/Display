START_GLOBALS
    DEF_GLOBAL( Alloc_checking_enabled, VIO_BOOL, FALSE )

    DEF_GLOBAL_COLOUR( Initial_background_colour )

    DEF_GLOBAL( Progress_threshold, int, 1000 )

    DEF_GLOBAL( Maximum_display_time, VIO_Real, 0.3 )
    DEF_GLOBAL( Min_interval_between_events, VIO_Real, 0.03 )
    DEF_GLOBAL( Min_interval_between_updates, VIO_Real, 0.02 )
    DEF_GLOBAL( Max_background_seconds, VIO_Real, 0.2 )
    DEF_GLOBAL( Event_timeout_factor, VIO_Real, 0.5 )
    DEF_GLOBAL( Size_of_interrupted, int, 250 )
    DEF_GLOBAL( Interval_of_check, int, 250 )

    DEF_GLOBAL( Initial_perspective_flag, VIO_BOOL, FALSE )
    DEF_GLOBAL( Perspective_distance_factor, VIO_Real, 2.0 )
    DEF_GLOBAL( Closest_front_plane, VIO_Real, 1.0e-5 )
    DEF_GLOBAL( Initial_render_mode, VIO_BOOL, TRUE )
    DEF_GLOBAL( Initial_shading_type, int, 1 )
    DEF_GLOBAL( Initial_light_switch, VIO_BOOL, TRUE )
    DEF_GLOBAL( Slice_double_buffer_flag, VIO_BOOL, TRUE )
    DEF_GLOBAL( Graphics_double_buffer_flag, VIO_BOOL, TRUE )
    DEF_GLOBAL( Graphics_transparency_flag, VIO_BOOL, TRUE )
    DEF_GLOBAL( Slice_readout_plane, int, 1 )
    DEF_GLOBAL( Display_frame_info, VIO_BOOL, FALSE )
    DEF_GLOBAL( Display_update_min, VIO_Real, 1.0 )
    DEF_GLOBAL( Frame_info_x, VIO_Real, 10.0 )
    DEF_GLOBAL( Frame_info_y, VIO_Real, 10.0 )
    DEF_GLOBAL( Initial_n_pixels_redraw, int, 100 )
    DEF_GLOBAL( Initial_slice_update_time, VIO_Real, 0.01 )
    DEF_GLOBAL( Initial_total_slice_update_time1, VIO_Real, 0.5 )
    DEF_GLOBAL( Initial_total_slice_update_time2, VIO_Real, 2.0 )
    DEF_GLOBAL( Initial_incremental_update, VIO_BOOL, FALSE )

    DEF_GLOBAL_COLOUR( Volume_outline_colour )
    DEF_GLOBAL_COLOUR( Cross_section_colour )
    DEF_GLOBAL5( Cross_section_spr, VIO_Surfprop, 0.3f, 0.6f, 0.6f, 40.0f, 1.0f )

    DEF_GLOBAL_COLOUR( Ambient_light_colour )
    DEF_GLOBAL_COLOUR( Directional_light_colour )
    DEF_GLOBAL3( Light_direction, VIO_Vector, 1.0f, -1.0f, -1.0f )

    DEF_GLOBAL( Initial_menu_window_width, int, 850 )
    DEF_GLOBAL( Initial_menu_window_height, int, 260 )
    DEF_GLOBAL( Canonical_menu_window_width, int, 850 )
    DEF_GLOBAL( Canonical_menu_window_height, int, 260 )
    DEF_GLOBAL( Menu_window_font, int, 1 )
    DEF_GLOBAL( Menu_window_font_size, VIO_Real, 10.0 )
    DEF_GLOBAL( X_menu_origin, VIO_Real, 70.0 )
    DEF_GLOBAL( X_menu_dx, VIO_Real, 10.0 )
    DEF_GLOBAL( X_menu_dy, VIO_Real, 0.0 )
    DEF_GLOBAL( X_menu_text_offset, VIO_Real, 5.0 )
    DEF_GLOBAL( Y_menu_origin, VIO_Real, 10.0 )
    DEF_GLOBAL( Y_menu_dx, VIO_Real, -15.0 )
    DEF_GLOBAL( Y_menu_dy, VIO_Real, 10.0 )
    DEF_GLOBAL( Y_menu_text_offset, VIO_Real, 15.0 )
    DEF_GLOBAL( Menu_n_chars_per_entry, int, 10 )
    DEF_GLOBAL( Menu_n_lines_per_entry, int, 2 )
    DEF_GLOBAL( Menu_character_height, VIO_Real, 20.0 )
    DEF_GLOBAL( Menu_character_width, VIO_Real, 7.5 )
    DEF_GLOBAL_COLOUR( Menu_character_colour )
    DEF_GLOBAL_COLOUR( Menu_character_inactive_colour )
    DEF_GLOBAL_COLOUR( Menu_box_colour )
    DEF_GLOBAL_COLOUR( Menu_key_colour )
    DEF_GLOBAL( Menu_key_character_offset, VIO_Real, 2.0 )
    DEF_GLOBAL( Initial_2_sided_flag, VIO_BOOL, FALSE )
    DEF_GLOBAL( Initial_backface_flag, VIO_BOOL, FALSE )
    DEF_GLOBAL( Initial_line_curves_flag, VIO_BOOL, FALSE )
    DEF_GLOBAL( Initial_n_curve_segments, int, 8 )

    DEF_GLOBAL( Segmenting_connectivity, int, 0 )

    DEF_GLOBAL( Visibility_on_input, VIO_BOOL, TRUE )

    DEF_GLOBAL( Initial_x_scale, VIO_Real, 1.0 )
    DEF_GLOBAL( Initial_y_scale, VIO_Real, 1.0 )
    DEF_GLOBAL( Initial_z_scale, VIO_Real, 1.0 )

    DEF_GLOBAL3( Default_line_of_sight, VIO_Vector, 0.0f, 0.0f, -1.0f )
    DEF_GLOBAL3( Default_horizontal, VIO_Vector, 1.0f, 0.0f, 0.0f )

    DEF_GLOBAL( Monitor_widths_to_eye, VIO_Real, 1.2 )

    DEF_GLOBAL( Default_marker_type, int, 0 )
    DEF_GLOBAL( Default_marker_structure_id, int, 1 )
    DEF_GLOBAL( Default_marker_patient_id, int, 1 )
    DEF_GLOBAL( Default_marker_size, VIO_Real, 1.0 )
    DEF_GLOBAL_COLOUR( Default_marker_colour )
    DEF_GLOBAL( Default_marker_label, VIO_STR, "Marker" )
    DEF_GLOBAL( Marker_pick_size, VIO_Real, 3.0 )

    DEF_GLOBAL_COLOUR( Viewport_feedback_colour )
    DEF_GLOBAL( Viewport_min_x_size, VIO_Real, 0.1 )
    DEF_GLOBAL( Viewport_min_y_size, VIO_Real, 0.1 )

    DEF_GLOBAL( Slice_magnification_step, VIO_Real, 1.2 )
    DEF_GLOBAL( Slice_fit_oversize, VIO_Real, 0.05 )
    DEF_GLOBAL_COLOUR( Slice_divider_colour )
    DEF_GLOBAL( Slice_divider_left, int, 5 )
    DEF_GLOBAL( Slice_divider_right, int, 5 )
    DEF_GLOBAL( Slice_divider_top, int, 5 )
    DEF_GLOBAL( Slice_divider_bottom, int, 5 )
    DEF_GLOBAL( Slice_divider_x_position, VIO_Real, 0.5 )
    DEF_GLOBAL( Slice_divider_y_position, VIO_Real, 0.5 )
    DEF_GLOBAL( Left_panel_width, int, 90 )
    DEF_GLOBAL( Text_panel_height, int, 190 )
    DEF_GLOBAL_COLOUR( Slice_text_colour )
    DEF_GLOBAL( Slice_text_font, int, 0 )
    DEF_GLOBAL( Slice_text_font_size, VIO_Real, 18.0 )
    DEF_GLOBAL( Slice_readout_text_font, int, 0 )
    DEF_GLOBAL( Slice_readout_text_font_size, VIO_Real, 12.0 )
    DEF_GLOBAL3( Slice_index_offset, VIO_Point, 10.0f, 5.0f, 0.0f )
    DEF_GLOBAL( Readout_text_colour, int, 2 )
    DEF_GLOBAL_COLOUR( Readout_text_rgb_colour )
    DEF_GLOBAL( Slice_index_x_format, VIO_STR, "X = %8g" )
    DEF_GLOBAL( Slice_index_y_format, VIO_STR, "Y = %8g" )
    DEF_GLOBAL( Slice_index_z_format, VIO_STR, "Z = %8g" )
    DEF_GLOBAL( Slice_probe_volume_index_format, VIO_STR, "V# %2d" )
    DEF_GLOBAL( Slice_probe_x_voxel_format, VIO_STR, "Xv %6.1f" )
    DEF_GLOBAL( Slice_probe_y_voxel_format, VIO_STR, "Yv %6.1f" )
    DEF_GLOBAL( Slice_probe_z_voxel_format, VIO_STR, "Zv %6.1f" )
    DEF_GLOBAL( Slice_probe_x_world_format, VIO_STR, "Xw %6.1f" )
    DEF_GLOBAL( Slice_probe_y_world_format, VIO_STR, "Yw %6.1f" )
    DEF_GLOBAL( Slice_probe_z_world_format, VIO_STR, "Zw %6.1f" )
    DEF_GLOBAL( Slice_probe_voxel_format,   VIO_STR, "Vx %6.5g" )
    DEF_GLOBAL( Slice_probe_val_format,     VIO_STR, "Vl %6.5g" )
    DEF_GLOBAL( Slice_probe_label_format,   VIO_STR, "Lb %6d" )
    DEF_GLOBAL( Slice_probe_distance_format, VIO_STR, "D%s %6.4g")
    DEF_GLOBAL( Probe_x_pos, int, 2 )
    DEF_GLOBAL( Probe_y_pos, int, 4 )
    DEF_GLOBAL( Probe_x_delta, int, 0 )
    DEF_GLOBAL( Probe_y_delta, int, 15 )
    DEF_GLOBAL( Slice_view1_axis1, int, 1 )
    DEF_GLOBAL( Slice_view1_axis2, int, 2 )

    DEF_GLOBAL( Slice_view2_axis1, int, 0 )
    DEF_GLOBAL( Slice_view2_axis2, int, 2 )

    DEF_GLOBAL( Slice_view3_axis1, int, 0 )
    DEF_GLOBAL( Slice_view3_axis2, int, 1 )

    DEF_GLOBAL3( Cursor_home, VIO_Point, 0.0f, 0.0f, 0.0f )
    DEF_GLOBAL_COLOUR( Cursor_colour )
    DEF_GLOBAL_COLOUR( Cursor_rgb_colour )
    DEF_GLOBAL( Show_cursor_contours, VIO_BOOL, FALSE )
    DEF_GLOBAL( Cursor_contour_overlay_flag, VIO_BOOL, FALSE )
    DEF_GLOBAL( Cursor_contour_thickness, VIO_Real, 3.0 )
    DEF_GLOBAL( Cursor_beep_on_surface, int, 0 )
    DEF_GLOBAL( Cursor_size_factor, VIO_Real, 0.05 )
    DEF_GLOBAL( Cursor_axis_size, VIO_Real, 50.0 )

    DEF_GLOBAL( Cursor_hor_start_0, VIO_Real, 3.0 )
    DEF_GLOBAL( Cursor_hor_end_0, VIO_Real, 10.0 )
    DEF_GLOBAL( Cursor_vert_start_0, VIO_Real, 3.0 )
    DEF_GLOBAL( Cursor_vert_end_0, VIO_Real, 10.0 )

    DEF_GLOBAL( Cursor_hor_start_1, VIO_Real, 3.0 )
    DEF_GLOBAL( Cursor_hor_end_1, VIO_Real, 10.0 )
    DEF_GLOBAL( Cursor_vert_start_1, VIO_Real, 3.0 )
    DEF_GLOBAL( Cursor_vert_end_1, VIO_Real, 10.0 )

    DEF_GLOBAL( Cursor_hor_start_2, VIO_Real, 3.0 )
    DEF_GLOBAL( Cursor_hor_end_2, VIO_Real, 30.0 )
    DEF_GLOBAL( Cursor_vert_start_2, VIO_Real, 3.0 )
    DEF_GLOBAL( Cursor_vert_end_2, VIO_Real, 30.0 )

    DEF_GLOBAL( Cursor_hor_start_3, VIO_Real, 10.0 )
    DEF_GLOBAL( Cursor_hor_end_3, VIO_Real, 30.0 )
    DEF_GLOBAL( Cursor_vert_start_3, VIO_Real, 10.0 )
    DEF_GLOBAL( Cursor_vert_end_3, VIO_Real, 30.0 )

    DEF_GLOBAL_COLOUR( Slice_cursor_colour1 )
    DEF_GLOBAL_COLOUR( Slice_cursor_colour2 )
    DEF_GLOBAL( Cursor_mouse_threshold, VIO_Real, 0.0 )
    DEF_GLOBAL( Max_cursor_angle, VIO_Real, 20.0 )
    DEF_GLOBAL( Min_cursor_angle_diff, VIO_Real, 5.0 )
    DEF_GLOBAL( Cursor_pick_distance, VIO_Real, 10.0 )

    DEF_GLOBAL_COLOUR( Slice_cross_section_colour )
    DEF_GLOBAL_COLOUR( Slice_crop_box_colour )
    DEF_GLOBAL( Slice_crop_pick_distance, VIO_Real, 10.0 )

    DEF_GLOBAL_COLOUR( Extracted_surface_colour )
    DEF_GLOBAL5( Default_surface_property, VIO_Surfprop, 0.3f, 0.3f, 0.4f, 10.0f, 1.0f )

    DEF_GLOBAL( Min_voxels_per_update, int, 10 )
    DEF_GLOBAL( Max_voxels_per_update, int, 100000 )
    DEF_GLOBAL( Max_seconds_per_voxel_update, VIO_Real, 0.1 )
    DEF_GLOBAL( Max_surface_refinements, int, 0 )
    DEF_GLOBAL( Max_surface_error, VIO_Real, 0.02 )
    DEF_GLOBAL( Edge_point_threshold, VIO_Real, 0.25 )
    DEF_GLOBAL( Edge_point_new_density, VIO_Real, 0.125 )
    DEF_GLOBAL( Voxel_validity_if_mixed, VIO_BOOL, FALSE )
    DEF_GLOBAL( Set_invalid_to_zero, VIO_BOOL, FALSE )

    DEF_GLOBAL( Max_fast_colour_lookup, int, 100000 )
    DEF_GLOBAL( Colour_table_size, int, 1000 )
    DEF_GLOBAL_COLOUR( Colour_below )
    DEF_GLOBAL( User_defined_n_intervals, int, 4 )
    DEF_GLOBAL_COLOUR( User_defined_min_colour )
    DEF_GLOBAL_COLOUR( User_defined_max_colour )
    DEF_GLOBAL( User_defined_colour_coding_flip, VIO_BOOL, TRUE )
    DEF_GLOBAL( Initial_low_absolute_position, VIO_Real, -1 )
    DEF_GLOBAL( Initial_high_absolute_position, VIO_Real, -1 )
    DEF_GLOBAL( Initial_histogram_contrast, VIO_BOOL, TRUE )
    DEF_GLOBAL( Initial_histogram_low_clip_index, int, 4 )
    DEF_GLOBAL( Initial_histogram_low, VIO_Real, 0.2 )
    DEF_GLOBAL( Initial_histogram_high, VIO_Real, 0.99999)
    DEF_GLOBAL( Initial_low_limit_position, VIO_Real, 0.25 )
    DEF_GLOBAL( Initial_high_limit_position, VIO_Real, 0.75 )
    DEF_GLOBAL_COLOUR( Colour_above )
    DEF_GLOBAL( Initial_colour_coding_type, int, 1 )
    DEF_GLOBAL( User_defined_interpolation_space, int, 1 )
    DEF_GLOBAL( Default_filter_type, int, 0 )
    DEF_GLOBAL( Default_filter_width, VIO_Real, 4.0 )

    DEF_GLOBAL( Colour_bar_resolution, int, 137 )
    DEF_GLOBAL( Colour_bar_top_offset, VIO_Real, 10.0 )
    DEF_GLOBAL( Colour_bar_bottom_offset, VIO_Real, 10.0 )
    DEF_GLOBAL( Colour_bar_left_offset, VIO_Real, 5.0 )
    DEF_GLOBAL( Colour_bar_width, VIO_Real, 10.0 )
    DEF_GLOBAL( Colour_bar_tick_width, VIO_Real, 3.0 )
    DEF_GLOBAL_COLOUR( Colour_bar_tick_colour )
    DEF_GLOBAL_COLOUR( Colour_bar_text_colour )
    DEF_GLOBAL( Colour_bar_text_font, int, 1 )
    DEF_GLOBAL( Colour_bar_text_size, VIO_Real, 10.0 )
    DEF_GLOBAL_COLOUR( Colour_bar_min_limit_colour )
    DEF_GLOBAL_COLOUR( Colour_bar_max_limit_colour )
    DEF_GLOBAL_COLOUR( Colour_bar_range_colour )
    DEF_GLOBAL( Colour_bar_desired_intervals, int, 5 )
    DEF_GLOBAL( Colour_bar_closest_text, VIO_Real, 10.0 )
    DEF_GLOBAL( Colour_bar_number_format, VIO_STR, "%g" )

    DEF_GLOBAL( Normal_towards_lower, VIO_BOOL, TRUE )

    DEF_GLOBAL( Initial_display_labels, VIO_BOOL, TRUE )
    DEF_GLOBAL( Initial_share_labels, VIO_BOOL, TRUE )
    DEF_GLOBAL( Initial_num_labels, int, 256 )
    DEF_GLOBAL( Label_colour_opacity, VIO_Real, 0.7 )

    DEF_GLOBAL( Use_cursor_origin, VIO_BOOL, TRUE )

    DEF_GLOBAL( One_active_flag, VIO_BOOL, FALSE )

    DEF_GLOBAL( Output_every, int, 0 )

    DEF_GLOBAL( Tmp_surface_name, VIO_STR, "" )

    DEF_GLOBAL( Marching_cubes_method, int, 1 )

    DEF_GLOBAL( Compute_neighbours_on_input, VIO_BOOL, FALSE )
    DEF_GLOBAL( Polygon_bintree_threshold, int, -1 )
    DEF_GLOBAL( Bintree_size_factor, VIO_Real, 0.3 )

    DEF_GLOBAL( Save_format, int, 0 )

    DEF_GLOBAL( Max_smoothing_distance, VIO_Real, 0.5 )
    DEF_GLOBAL( Smoothing_ratio, VIO_Real, 0.5 )
    DEF_GLOBAL( Smoothing_normal_ratio, VIO_Real, 0.0 )
    DEF_GLOBAL( Smoothing_threshold, VIO_Real, 0.01 )

    DEF_GLOBAL( Cursor_bitplanes, int, 0 )

    DEF_GLOBAL_COLOUR( Surface_curve_colour )
    DEF_GLOBAL( Surface_curve_overlay_flag, VIO_BOOL, FALSE )
    DEF_GLOBAL( Surface_curve_thickness, VIO_Real, 3.0 )
    DEF_GLOBAL( Line_curvature_weight, VIO_Real, 0.0 )
    DEF_GLOBAL( Min_surface_curve_curvature, VIO_Real, 0.0 )
    DEF_GLOBAL( Max_surface_curve_curvature, VIO_Real, 0.0 )

    DEF_GLOBAL( N_painting_polygons, int, 100 )
    DEF_GLOBAL_COLOUR( Visible_segmenting_colour )
    DEF_GLOBAL_COLOUR( Invisible_segmenting_colour )
    DEF_GLOBAL( Snap_to_polygon_vertex, VIO_BOOL, TRUE )

    DEF_GLOBAL( N_fitting_samples, int, 2000 )
    DEF_GLOBAL( Isovalue_factor, VIO_Real, 0.0 )
    DEF_GLOBAL( Fitting_isovalue, VIO_Real, 1000.0 )
    DEF_GLOBAL( Gradient_strength_factor, VIO_Real, 1.0 )
    DEF_GLOBAL( Gradient_strength_exponent, VIO_Real, 1.0 )
    DEF_GLOBAL( Curvature_factor, VIO_Real, 0.5 )
    DEF_GLOBAL( Surface_point_distance_threshold, VIO_Real, 1.0 )
    DEF_GLOBAL( Surface_point_distance_factor, VIO_Real, 1.0 )
    DEF_GLOBAL( Fitting_tolerance, VIO_Real, 1.0e-4 )
    DEF_GLOBAL( Max_fitting_evaluations, int, 50 )
    DEF_GLOBAL( Surface_model_resolution, int, 30 )
    DEF_GLOBAL( Minimization_method, int, 1 )
    DEF_GLOBAL( Max_parameter_delta, VIO_Real, 1.0 )

    DEF_GLOBAL( Max_voxel_scan_distance, VIO_Real, 2.0 )
    DEF_GLOBAL( Max_parametric_scan_distance, VIO_Real, 0.1 )
    DEF_GLOBAL( Max_polygon_scan_distance, VIO_Real, 2.0 )


    DEF_GLOBAL( Initial_atlas_opacity, VIO_Real, 1.0 )
    DEF_GLOBAL( Initial_atlas_transparent_threshold, int, 220 )
    DEF_GLOBAL( Initial_atlas_tolerance_x, VIO_Real, 1.0 )
    DEF_GLOBAL( Initial_atlas_tolerance_y, VIO_Real, 1.0 )
    DEF_GLOBAL( Initial_atlas_tolerance_z, VIO_Real, 1.0 )
    DEF_GLOBAL( Atlas_filename, VIO_STR, "/avgbrain/atlas/talairach/obj/Talairach_atlas.list" )
    DEF_GLOBAL( Default_atlas_state, VIO_BOOL, FALSE )

    DEF_GLOBAL( Default_x_brush_radius, VIO_Real, 3.0 )
    DEF_GLOBAL( Default_y_brush_radius, VIO_Real, 3.0 )
    DEF_GLOBAL( Default_z_brush_radius, VIO_Real, 0.0 )
    DEF_GLOBAL( Default_paint_label, int, 1 )
    DEF_GLOBAL( Default_fast_painting_flag, VIO_BOOL, TRUE )
    DEF_GLOBAL( Default_cursor_follows_paintbrush_flag, VIO_BOOL, FALSE )
    DEF_GLOBAL( Initial_mouse_scale_factor, VIO_Real, 1.0 )
    DEF_GLOBAL( Draw_brush_outline, VIO_BOOL, TRUE )
    DEF_GLOBAL_COLOUR( Brush_outline_colour )
    DEF_GLOBAL( Brush_outline_offset, int, 0 )
    DEF_GLOBAL( Snap_brush_to_centres, VIO_BOOL, TRUE )
    DEF_GLOBAL( Undo_enabled, VIO_BOOL, TRUE )

    DEF_GLOBAL( Volume_continuity, int, 0 )
    DEF_GLOBAL( Convert_volumes_to_byte, VIO_BOOL, TRUE )

    DEF_GLOBAL( Move_slice_speed, VIO_Real, 0.25 )
    DEF_GLOBAL( Pixels_per_double_size, VIO_Real, 100.0 )

    DEF_GLOBAL( Clear_before_polygon_scan, VIO_BOOL, FALSE )

    DEF_GLOBAL_COLOUR( Histogram_colour )
    DEF_GLOBAL( Histogram_x_scale, VIO_Real, 1.0 )
    DEF_GLOBAL( Histogram_smoothness_ratio, VIO_Real, 0.02 )
    DEF_GLOBAL( Histogram_extra_width, VIO_Real, 0.0 )

    DEF_GLOBAL_COLOUR( Menu_name_colour )
    DEF_GLOBAL( Menu_name_x, VIO_Real, 5.0 )
    DEF_GLOBAL( Menu_name_y, VIO_Real, 20.0 )
    DEF_GLOBAL( Menu_name_font, int, 1 )
    DEF_GLOBAL( Menu_name_font_size, VIO_Real, 12.0 )

    DEF_GLOBAL( Use_transparency_hardware, VIO_BOOL, TRUE )

    DEF_GLOBAL( Crop_volume_command, VIO_STR, "mincreshape %s %s -start %d,%d,%d -count %d,%d,%d" )

    DEF_GLOBAL( Crop_volumes_on_input, VIO_BOOL, FALSE )
    DEF_GLOBAL( Crop_if_smaller, VIO_Real, 0.75 )

    DEF_GLOBAL( Initial_slice_continuity, int, -1 )

    DEF_GLOBAL_COLOUR( Cursor_pos_colour )
    DEF_GLOBAL( Cursor_pos_x_origin, VIO_Real, 400.0 )
    DEF_GLOBAL( Cursor_pos_y_origin, VIO_Real, 30.0 )
    DEF_GLOBAL( Cursor_pos_title, VIO_STR, "Cursor Pos:" )
    DEF_GLOBAL( Cursor_pos_format, VIO_STR, "%.1f %.1f %.1f" )

    DEF_GLOBAL_COLOUR( Unfinished_flag_colour )
    DEF_GLOBAL( Unfinished_flag_width, VIO_Real, 5.0 )

    DEF_GLOBAL( Initial_volume_rotation_step, VIO_Real, 1.0 )
    DEF_GLOBAL( Initial_volume_translation_step, VIO_Real, 1.0 )
    DEF_GLOBAL( Initial_volume_scale_step, VIO_Real, 1.05 )

    DEF_GLOBAL( Enable_volume_caching, VIO_BOOL, TRUE )
    DEF_GLOBAL( Volume_cache_threshold, int, -1 )
    DEF_GLOBAL( Volume_cache_size, int, -1 )
    DEF_GLOBAL( Volume_cache_block_size, int, -1 )

    DEF_GLOBAL( Crop_label_volumes_threshold, VIO_Real, 0.9 )
    DEF_GLOBAL( Initial_crop_labels_on_output, VIO_BOOL, FALSE )

    DEF_GLOBAL( Duplicate_boundary_vertices, VIO_BOOL, FALSE )

    DEF_GLOBAL( Toggle_freestyle_painting, VIO_BOOL, TRUE )
    DEF_GLOBAL( First_straightline_right_mouse_down,VIO_BOOL, TRUE )

    DEF_GLOBAL( Output_label_filename, VIO_STR, "" )
    DEF_GLOBAL( Exit_error_load_file, VIO_BOOL, TRUE )
    DEF_GLOBAL( Ratio_volume_index, VIO_STR, "" )
    DEF_GLOBAL( Ratio_volume_index_format, VIO_STR, "%d,%d" )
    DEF_GLOBAL( Slice_probe_ratio_format, VIO_STR, "Ratio %6.2f" )
    DEF_GLOBAL_COLOUR( Slice_probe_ratio_colour )

    DEF_GLOBAL( Hide_3D_window, VIO_BOOL, TRUE)
    DEF_GLOBAL( Hide_marker_window, VIO_BOOL, TRUE)
    DEF_GLOBAL( Hide_menu_window, VIO_BOOL, FALSE)

    DEF_GLOBAL( Initial_marker_window_width, int, 250 )
    DEF_GLOBAL( Initial_marker_window_height, int, 500 )
    DEF_GLOBAL( Canonical_marker_window_width, int, 250 )
    DEF_GLOBAL( Canonical_marker_window_height, int, 500 )

    DEF_GLOBAL( N_selected_displayed, int, 50 )
    DEF_GLOBAL( Selected_x_origin, VIO_Real, 20.0 )
    DEF_GLOBAL( Selected_y_origin, VIO_Real, 480.0 )
    DEF_GLOBAL_COLOUR( Selected_colour )
    DEF_GLOBAL_COLOUR( Visible_colour )
    DEF_GLOBAL_COLOUR( Invisible_colour )
    DEF_GLOBAL( Character_height_in_pixels, VIO_Real, 10.0 )
    DEF_GLOBAL( Selected_box_x_offset, VIO_Real, 2.0 )
    DEF_GLOBAL( Selected_box_y_offset, VIO_Real, 2.0 )
    DEF_GLOBAL( Tags_from_label, VIO_BOOL, FALSE )

    DEF_GLOBAL( Initial_undo_feature, VIO_BOOL, TRUE )
  
    DEF_GLOBAL( Object_outline_width, VIO_Real, 1.0 )
    DEF_GLOBAL( Object_outline_enabled, VIO_BOOL, TRUE )

END_GLOBALS
