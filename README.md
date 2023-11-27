# Submodel
A command-line tool to automatically derive submodels from a threedi model.

## Features
Some features are:
* Automatically derived for a set of subarea (polygon) geometries
* Clips all 1D elements
* Clips all 2D rasters
* Keep alignment of grid cells due to raster outlining with calculation grid cells of original model
* Verifies input for tool before running

## Input
Required:
* Threedi schematisation directory; folder including .sqlite file, .gpkg file and rasters folder.
* Subareas; path to vector containing geometry of submodel areas as different features. Most common shape file extensions supported.
* Field name; name of field containing names of submodels. This name will be used to name the output folder and the submodel. Should contain unique subarea names.
* Calculation grid cells; path to vector containing geometry of calculation grid cells. Most common shape file extensions supported.

Optional:
* Subareas layer name; will be used as layer name in case subareas_path is of format geopackage (.gpkg).
* Calculation grid cells layer name; will be used as layer name in case calculation_grid_cells_path is of format geopackage (.gpkg).

## Usage
The step-by-step approach is:
1. Load the original schematisation (sqlite) in the Schematisation Editor (QGIS / MI). This will generate a 3Di schematisation geopackage.
2. Use the processing tool of the Schematisation Editor to generate the grid from schematisation. This will generate the calculation grid cells.
3. Create or open the subarea polygon vector layer. Open the attribute table. Make sure there is a field with subarea names. Make sure the names are unique.
4. Run the tool.

```
submodel.py [-h] [--subareas_layer_name SUBAREAS_LAYER_NAME]
                   [--calculation_grid_cells_layer_name CALCULATION_GRID_CELLS_LAYER_NAME]
                   schematisation_directory subareas_path field_name calculation_grid_cells_path

positional arguments:
  schematisation_directory
                        Required. Folder including .sqlite file, .gpkg file and rasters folder.
  subareas_path         Required. Path to vector containing geometry of submodel areas as different features. Most
                        common shape file extensions supported.
  field_name            Required. Name of field containing names of submodels. This name will be used to name output
                        folder and submodel. Should contain unique subarea names.
  calculation_grid_cells_path
                        Required. Path to vector containing geometry of calculation grid cells. Most common shape file
                        extensions supported.

options:
  -h, --help            show this help message and exit
  --subareas_layer_name SUBAREAS_LAYER_NAME, -s SUBAREAS_LAYER_NAME
                        Optional. Will be used as layer name in case subareas_path is of format geopackage (.gpkg)
  --calculation_grid_cells_layer_name CALCULATION_GRID_CELLS_LAYER_NAME, -c CALCULATION_GRID_CELLS_LAYER_NAME
                        Optional. Will be used as layer name in case calculation_grid_cells_path is of format
                        geopackage (.gpkg)
```

5. Open the 3Di schematisation geopackage in the Schematisation Editor.
6. Save to spatiallite. 

Now you are ready to upload the schematisation to the 3Di cloud.

## Dependencies
[Dependencies](dep)



