/**
 * \file poly_formats.c
 * \brief Functions to write polygonal objects in alternate formats.
 *
 * Note that these functions have been intentionally kept simple - we
 * want to write files that can be readily imported into tools such as
 * Blender, but we do _not_ care about handling every conceivable
 * case.
 *
 * One important omission is standardizing the coordinate space where
 * possible. It is not clear what the standard coordinate space is for
 * most of these formats. X3D probably defines something along these
 * lines, but I don't implement it here (yet).
 */

#include <display.h>

#if GIFTI_FOUND
#include <time.h>
#include "gifti_io.h"
#endif /* GIFTI_FOUND */
/**
 * Write a Wavefront .obj file in ASCII format, given a memory representation
 * of an MNI polygonal surface.
 *
 * This function will write a SINGLE polygonal surface object, with 
 * o colour information. These files do not support vertex or face colouring.
 *
 * References:
 * https://en.wikipedia.org/wiki/Wavefront_.obj_file
 *
 * \param filename The name of the file to create.
 * \param object_ptr A pointer to a single polygonal object.
 * \returns VIO_OK if the operation succeeded.
 */
VIO_Status
output_wavefront_obj(VIO_STR filename, object_struct *object_ptr)
{
    polygons_struct *polygons_ptr = get_polygons_ptr( object_ptr );
    FILE *fp = fopen(filename, "w");
    int i = 0, j = 0;

    if ( fp == NULL || polygons_ptr == NULL )
    {
        return VIO_ERROR;
    }

    fprintf( fp, "# Created by %s %s\n", PROJECT_NAME, PROJECT_VERSION );
    for (i = 0; i < polygons_ptr->n_points; i++)
    {
        fprintf( fp, "v %f %f %f\n",
                 Point_x(polygons_ptr->points[i]),
                 Point_y(polygons_ptr->points[i]),
                 Point_z(polygons_ptr->points[i]) );
    }
    if ( polygons_ptr->normals != NULL )
    {
        for (i = 0; i < polygons_ptr->n_points; i++)
        {
            fprintf( fp, "vn %f %f %f\n",
                     Vector_x(polygons_ptr->normals[i]),
                     Vector_y(polygons_ptr->normals[i]),
                     Vector_z(polygons_ptr->normals[i]) );
        }
    }
    
    for (j = i = 0; i < polygons_ptr->n_items; i++)
    {
        int k = polygons_ptr->end_indices[i];
        fprintf( fp, "f " );
        for ( ; j < k; j++)
        {
            int n = polygons_ptr->indices[j] + 1;
            fprintf( fp, "%d//%d", n, n );
            if (j < k - 1)
                fprintf( fp, " " );
        }
        fprintf( fp, "\n");
    }
    fclose( fp );
    return VIO_OK;
}

/**
 * Write a Stanford .ply file in ASCII format, given a memory representation
 * of an MNI polygonal surface.
 *
 * This function will write a SINGLE polygonal surface object, with either
 * a no colour information, or with per-vertex colours. It isn't trivial for
 * the Stanford format to represent a per-face colour.
 *
 * References:
 * https://en.wikipedia.org/wiki/PLY_(file_format)
 * http://paulbourke.net/dataformats/ply/
 *
 * \param filename The name of the file to create.
 * \param object_ptr A pointer to a single polygonal object.
 * \returns VIO_OK if the operation succeeded.
 */
VIO_Status
output_stanford_ply(VIO_STR filename, object_struct *object_ptr)
{
    polygons_struct *polygons_ptr = get_polygons_ptr( object_ptr );
    FILE *fp = fopen(filename, "w");
    int i = 0, j = 0;

    if ( fp == NULL || polygons_ptr == NULL )
    {
        return VIO_ERROR;
    }

    fprintf( fp, "ply\n");
    fprintf( fp, "format ascii 1.0\n");
    fprintf( fp, "comment Created by %s %s\n", PROJECT_NAME, PROJECT_VERSION );
    fprintf( fp, "element vertex %d\n", polygons_ptr->n_points );
    fprintf( fp, "property float x\n" );
    fprintf( fp, "property float y\n" );
    fprintf( fp, "property float z\n" );
    if (polygons_ptr->colour_flag == PER_VERTEX_COLOURS)
    {
      fprintf( fp, "property uchar red\n" );
      fprintf( fp, "property uchar green\n" );
      fprintf( fp, "property uchar blue\n" );
    }
    fprintf( fp, "element face %d\n", polygons_ptr->n_items );
    fprintf( fp, "property list uchar int vertex_index\n" );
    fprintf( fp, "end_header\n" );
    if (polygons_ptr->colour_flag == PER_VERTEX_COLOURS)
    {
        for (i = 0; i < polygons_ptr->n_points; i++)
        {
            fprintf(fp, "%f %f %f %d %d %d\n",
                    Point_x(polygons_ptr->points[i]),
                    Point_y(polygons_ptr->points[i]),
                    Point_z(polygons_ptr->points[i]),
                    get_Colour_r(polygons_ptr->colours[i]),
                    get_Colour_g(polygons_ptr->colours[i]),
                    get_Colour_b(polygons_ptr->colours[i]));
        }
    }
    else
    {
        for (i = 0; i < polygons_ptr->n_points; i++)
        {
            fprintf(fp, "%f %f %f\n",
                    Point_x(polygons_ptr->points[i]),
                    Point_y(polygons_ptr->points[i]),
                    Point_z(polygons_ptr->points[i]));
        }
    }

    
    for (j = i = 0; i < polygons_ptr->n_items; i++)
    {
        int k = polygons_ptr->end_indices[i];
        fprintf( fp, "%d ", k - j );
        for ( ; j < k; j++)
        {
            fprintf( fp, "%d", polygons_ptr->indices[j] );
            if (j < k - 1)
                fprintf( fp, " " );
        }
        fprintf( fp, "\n");
    }
    fclose( fp );
    return VIO_OK;
}

/**
 * Write an X3D format file, given a memory representation
 * of an MNI polygonal surface.
 *
 * As with most XML formats, reading the actual specification makes my
 * head hurt. I relied on reverse-engineering to get this compatible
 * with Blender's import function.
 *
 * References:
 * https://en.wikipedia.org/wiki/X3D
 *
 * \param filename The name of the file to create.
 * \param object_ptr A pointer to a single polygonal object.
 * \returns VIO_OK if the operation succeeded.
 */
VIO_Status
output_x3d(VIO_STR filename, object_struct *object_ptr)
{
    polygons_struct *polygons_ptr = get_polygons_ptr( object_ptr );
    FILE *fp = fopen( filename, "w" );
    int i = 0, j = 0;

    if ( fp == NULL || polygons_ptr == NULL )
    {
        return VIO_ERROR;
    }

    fprintf( fp, "<?xml version='1.0' encoding='UTF-8'?>\n" );
    fprintf( fp, "<!DOCTYPE X3D PUBLIC 'ISO//Web3D//DTD X3D 3.0//EN' 'http://www.web3d.org/specifications/x3d-3.0.dtd'>\n");
    fprintf( fp, "<!-- Created by %s %s -->\n", PROJECT_NAME, PROJECT_VERSION );
    fprintf( fp, "<X3D><Scene><Transform><Shape>\n" );

    if (polygons_ptr->colour_flag == PER_VERTEX_COLOURS)
    {
        fprintf( fp, "<IndexedFaceSet colorPerVertex='true' ");
    }
    else
    {
        fprintf( fp, "<IndexedFaceSet colorPerVertex='false' ");
    }
    fprintf( fp, "coordIndex='");
    for (j = i = 0; i < polygons_ptr->n_items; i++)
    {
        int k = polygons_ptr->end_indices[i];
        for ( ; j < k; j++)
        {
            fprintf( fp, "%d ", polygons_ptr->indices[j] );
        }
        fprintf( fp, "-1");
        if (i < polygons_ptr->n_items - 1)
            fprintf( fp, " " );
    }
    if (polygons_ptr->colour_flag == PER_VERTEX_COLOURS)
    {
        fprintf( fp, "' colorIndex='" );
        for ( i = 0; i < polygons_ptr->n_points; i++ )
        {
            fprintf( fp, "%d", i );
            if (i < polygons_ptr->n_points - 1)
                fprintf( fp, " " );
        }
    }
    fprintf( fp, "'>\n" );

    fprintf( fp, "<Coordinate point='" );
    for (i = 0; i < polygons_ptr->n_points; i++)
    {
        fprintf(fp, "%f %f %f",
                Point_x(polygons_ptr->points[i]),
                Point_y(polygons_ptr->points[i]),
                Point_z(polygons_ptr->points[i]));
        if ( i < polygons_ptr->n_points - 1 )
            fprintf( fp, " " );
    }
    fprintf( fp, "'/>\n");

    if (polygons_ptr->colour_flag == PER_VERTEX_COLOURS )
    {
        fprintf( fp, "<Color color='" );
        for (i = 0; i < polygons_ptr->n_points; i++)
        {
            fprintf(fp, "%f %f %f",
                    get_Colour_r_0_1(polygons_ptr->colours[i]),
                    get_Colour_g_0_1(polygons_ptr->colours[i]),
                    get_Colour_b_0_1(polygons_ptr->colours[i]) );

            if ( i < polygons_ptr->n_points - 1 )
              fprintf( fp, " " );
        }
        fprintf( fp, "'/>");
    }
    fprintf( fp, "</IndexedFaceSet>\n");
    fprintf( fp, "</Shape></Transform></Scene></X3D>\n");
    fclose( fp );
    return VIO_OK;
}

#if GIFTI_FOUND
/**
 * Write a GIFTI format file for a surface. The GIfTI format allows for
 * the creation of a wide range of possible structures, but we are only
 * interested in the "Surface" file, which represents a triangular mesh.
 * This type of file does not represent colours or vertex values, although
 * other GIFTI file types can separately represent those things.
 *
 * \param filename The name of the file to create.
 * \param object_ptr A pointer to a single polygonal object.
 * \returns VIO_OK if the operation succeeded.
 */
VIO_Status
output_gifti(VIO_STR filename, object_struct *object_ptr)
{
    polygons_struct *polygons_ptr = get_polygons_ptr( object_ptr );
    gifti_image *gifti_ptr;
    int n;
    int i;
    float *point_data;
    int *index_data;
    time_t t;
    char time_str[128];

    /* Check that this is in fact a triangular mesh. Give up and 
     * indicate failure if not.
     */
    n = 0;
    for (i = 0; i < polygons_ptr->n_items; i++)
    {
        if (polygons_ptr->end_indices[i] - n != 3)
        {
            return VIO_ERROR;
        }
        n = polygons_ptr->end_indices[i];
    }

    /* Create the skeleton of the GIFTI image. This does not fill in any
     * of the details. We have to do that ourselves in the current 
     * version of the API.
     */
    gifti_ptr = gifti_create_image( 0, 0, 0, 0, NULL, 0);
    if (gifti_ptr == NULL)
    {
        return VIO_ERROR;
    }

    /* Add some potentially useful metadata to the file.
     */
    if ( gifti_add_to_meta( &gifti_ptr->meta, "MNI-Display-Version",
                            PACKAGE_VERSION , 0 ))
    {
        return VIO_ERROR;
    }

    t = time(NULL);
    strftime(time_str, sizeof(time_str) - 1, "%F %T", localtime(&t));
    if ( gifti_add_to_meta( &gifti_ptr->meta, "Date", time_str, 0))
    {
        return VIO_ERROR;
    }

    /* Now start setting up the data arrays that will represent the
     * actual surface data. The surface file format requires two 
     * data arrays, one for the points and another for the triangle
     * indices. They have to be in the right order.
     */
    if ( gifti_add_empty_darray(gifti_ptr, 2) )
    {
        return VIO_ERROR;
    }

    /* Fill in the data array for the points. */
    n = 0;
    gifti_set_DA_defaults(gifti_ptr->darray[n]);
    gifti_set_atr_in_DAs(gifti_ptr, "Intent", "NIFTI_INTENT_POINTSET", &n, 1);
    gifti_set_atr_in_DAs(gifti_ptr, "DataType", "NIFTI_TYPE_FLOAT32", &n, 1);
    gifti_ptr->darray[n]->num_dim = 2;
    gifti_ptr->darray[n]->dims[0] = polygons_ptr->n_points;
    gifti_ptr->darray[n]->dims[1] = 3;

    /* Fill in the data array for the indices. */
    n = 1;
    gifti_set_DA_defaults(gifti_ptr->darray[n]);
    gifti_set_atr_in_DAs(gifti_ptr, "Intent", "NIFTI_INTENT_TRIANGLE", &n, 1);
    gifti_set_atr_in_DAs(gifti_ptr, "DataType", "NIFTI_TYPE_INT32", &n, 1);
    gifti_ptr->darray[n]->num_dim = 2;
    gifti_ptr->darray[n]->dims[0] = polygons_ptr->n_items;
    gifti_ptr->darray[n]->dims[1] = 3;

    /* Now we need to update the number of bytes per item in the 
     * structure.
     */
    if (gifti_update_nbyper(gifti_ptr))
    {
        return VIO_ERROR;
    }

    /* And then update the number of values in the structure.
     */
    gifti_ptr->darray[0]->nvals = gifti_darray_nvals(gifti_ptr->darray[0]);
    gifti_ptr->darray[1]->nvals = gifti_darray_nvals(gifti_ptr->darray[1]);

    /* Now allocate the data we will need to store the actual data points.
     */
    if (gifti_alloc_DA_data(gifti_ptr, NULL, 0))
    {
        return VIO_ERROR;
    }

    /* Copy the points to the first data array. */
    point_data = (float *) gifti_ptr->darray[0]->data;
    for_less( i, 0, polygons_ptr->n_points )
    {
        point_data[i * 3 + 0] = Point_x(polygons_ptr->points[i]);
        point_data[i * 3 + 1] = Point_y(polygons_ptr->points[i]);
        point_data[i * 3 + 2] = Point_z(polygons_ptr->points[i]);
    }

    /* Copy the indices to the second data array. */
    index_data = (int *) gifti_ptr->darray[1]->data;
    for_less( i, 0, polygons_ptr->n_items * 3)
    {
        index_data[i] = polygons_ptr->indices[i];
    }

    /* Actually write the file. */
    gifti_write_image(gifti_ptr, filename, 1);

    /* Finally free the memory we allocated. */
    gifti_free_image(gifti_ptr);

    return VIO_OK;
}

#endif /* GIFTI_FOUND */
