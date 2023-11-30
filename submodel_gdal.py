import argparse
from pathlib import Path
import shutil
from tqdm import tqdm

import geopandas as gpd
from shapely.geometry import mapping
from shapely import unary_union
import fiona

from osgeo import gdal, osr
import json
import tempfile

import warnings
warnings.filterwarnings("ignore", category=UserWarning, message="You are attempting to write an empty DataFrame to file. For some drivers, this operation may fail.")

### Custom exceptions

class GeoPackageFileNotFoundError(Exception):
    pass

class SQLiteFileNotFoundError(Exception):
    pass

class RastersDirectoryNotFoundError(Exception):
    pass

class FileNotFoundError(Exception):
    pass

class LayerNotFoundError(Exception):
    pass

class SubareaLayerEmpty(Exception):
    pass

class FieldNameNotFoundError(Exception):
    pass

class SubareaNamesNotUnique(Exception):
    pass

class NoCalcGridCellsSelected(Exception):
    pass


class Submodels:
    name = "Threedi SubModel"
    
    def __init__(self, schematisation_directory, subareas_path, field_name, calculation_grid_cells_path, 
                 subareas_layer_name=None, calculation_grid_cells_layer_name=None):
        
        self.schematisation_directory = Path(schematisation_directory)
        self.subareas_path = Path(subareas_path)
        self.subareas_layer_name = subareas_layer_name
        self.field_name = field_name
        self.calculation_grid_cells_path = Path(calculation_grid_cells_path)
        self.calculation_grid_cells_layer_name = calculation_grid_cells_layer_name
        
        # Find files and check properties
        self.schematisation_gpkg = self.find_gpkg_file()
        self.schematisation_sqlite = self.find_sqlite_file()
        self.rasters_directory = self.find_rasters_directory()
        
        self.check_file_exists(self.subareas_path)
        self.subareas = self.read_to_gdf(self.subareas_path, self.subareas_layer_name)
        self.check_field_existence_and_uniqueness()
        
        self.check_file_exists(self.calculation_grid_cells_path)
        self.calculation_grid_cells = self.read_to_gdf(self.calculation_grid_cells_path, self.calculation_grid_cells_layer_name)
        
        # Clip
        for _, subarea in tqdm(self.subareas.iterrows(), total=len(self.subareas),
                                  desc="Clipping", unit="subarea"):
            self.clip(subarea)
        
    
    def find_gpkg_file(self):
        gpkg_files = self.schematisation_directory.glob('*.gpkg')
        if gpkg_files:
            return next(gpkg_files)
        else:
            raise GeoPackageFileNotFoundError("No GeoPackage file found in the schematisation directory.")    

    def find_sqlite_file(self):
        sqlite_files = self.schematisation_directory.glob('*.sqlite')
        if sqlite_files:
            return next(sqlite_files)
        else:
            raise SQLiteFileNotFoundError("No SQLite file found in the schematisation directory.")   

    def find_rasters_directory(self):
        for item in self.schematisation_directory.iterdir():
            if item.is_dir() and item.name == "rasters":
                return item
        
        raise RastersDirectoryNotFoundError("No folder 'rasters' found in the schematisation directory.")

    def check_file_exists(self, file_path):
        if file_path.exists():
            return 
        else:
            raise FileNotFoundError(f"The file at '{file_path}' does not exist.")

    def read_to_gdf(self, path, layer_name):
        if path.suffix == '.gpkg':
            try:
                gdf = gpd.read_file(path, layer=layer_name)
            except:
                raise LayerNotFoundError(f"The file at '{path}' does not contain a layer named '{layer_name}'.")
        else:
            gdf = gpd.read_file(path)
            
        return gdf
            
    def check_field_existence_and_uniqueness(self):
        if self.field_name in self.subareas.columns:
            if self.subareas.empty:
                raise SubareaLayerEmpty(f"The file at '{self.subareas_path}' does not contain any features.")
            if self.subareas[self.field_name].duplicated().any():
                raise SubareaNamesNotUnique(f"The values in the '{self.field_name}' column are not unique.")
        else:
            raise FieldNameNotFoundError(f"The file at '{self.subareas_path}' does not contain a field with the name '{self.field_name}'.")


    def clip(self, subarea):
        
        # Create a folder for the clipped schematisation
        name = subarea[self.field_name]
        output_directory = self.schematisation_directory / name
        output_directory.mkdir(parents=True, exist_ok=True)
        output_schematisation_gpkg_path = output_directory / (self.schematisation_gpkg.stem + '_' + name + self.schematisation_gpkg.suffix)
        output_schematisation_sqlite_path = output_directory / (self.schematisation_sqlite.stem + '_' + name + self.schematisation_sqlite.suffix)

        # Copy the schematisation to the new folder
        shutil.copy(self.schematisation_gpkg, output_schematisation_gpkg_path)
        shutil.copy(self.schematisation_sqlite, output_schematisation_sqlite_path)
        
        ### vector clipping
        # Read schematisation
        layers_dict = self.read_geopackage_layers(output_schematisation_gpkg_path)
        
        # Read clipping layers
        connection_node = layers_dict.get('connection_node')
        manhole = layers_dict.get('manhole')
        pipe = layers_dict.get('pipe')
        weir = layers_dict.get('weir')
        orifice = layers_dict.get('orifice')
        culvert = layers_dict.get('culvert')
        cross_section_location = layers_dict.get('cross_section_location')
        channel = layers_dict.get('channel')
        pumpstation_map = layers_dict.get('pumpstation_map')
        pumpstation = layers_dict.get('pumpstation')
        boundary_condition_1d = layers_dict.get('1d_boundary_condition')    #
        boundary_condition_2d = layers_dict.get('2d_boundary_condition')    #
        lateral_1d = layers_dict.get('1d_lateral')                          #
        lateral_2d = layers_dict.get('2d_lateral')                          #
        impervious_surface_map = layers_dict.get('impervious_surface_map')
        impervious_surface = layers_dict.get('impervious_surface')
        linear_obstacle = layers_dict.get('linear_obstacle')
        potential_breach = layers_dict.get('potential_breach')
        exchange_line = layers_dict.get('exchange_line')
        grid_refinement = layers_dict.get('grid_refinement')
        grid_refinement_area = layers_dict.get('grid_refinement_area')
        
        # subarea from series to GeoDataFrame
        subarea_gdf = gpd.GeoDataFrame(subarea.to_frame().T, geometry='geometry')
        subarea_gdf.crs = self.subareas.crs

        # Perform a spatial join between connection nodes and subarea
        filtered_connection_node = self.spatial_join(connection_node, subarea_gdf, 'inner', 'intersects', '_subarea')
        
        # Extract the IDs from the result
        valid_connection_node_ids = filtered_connection_node['id']
        
        # A conncetion node should always be connected to manhole, channel, pipe, orifice, culvert, weir, pumpstation
        
        # Filter 1d point items on a single connection node
        filtered_manhole = manhole[
            manhole['connection_node_id'].isin(valid_connection_node_ids)
        ]
        
        filtered_pumpstation = pumpstation[
            pumpstation['connection_node_id'].isin(valid_connection_node_ids)
        ]
        
        filtered_pipe = pipe[
            (pipe['connection_node_start_id'].isin(valid_connection_node_ids)) &
            (pipe['connection_node_end_id'].isin(valid_connection_node_ids))
        ]
        
        filtered_weir = weir[
            (weir['connection_node_start_id'].isin(valid_connection_node_ids)) &
            (weir['connection_node_end_id'].isin(valid_connection_node_ids))
        ]
        
        filtered_orifice = orifice[
            (orifice['connection_node_start_id'].isin(valid_connection_node_ids)) &
            (orifice['connection_node_end_id'].isin(valid_connection_node_ids))
        ]
        
        filtered_culvert = culvert[
            (culvert['connection_node_start_id'].isin(valid_connection_node_ids)) &
            (culvert['connection_node_end_id'].isin(valid_connection_node_ids))
        ]
        
        filtered_pumpstation_map = pumpstation_map[
            (pumpstation_map['connection_node_start_id'].isin(valid_connection_node_ids)) &
            (pumpstation_map['connection_node_end_id'].isin(valid_connection_node_ids))
        ]
        
        filtered_channel = channel[
            (channel['connection_node_start_id'].isin(valid_connection_node_ids)) &
            (channel['connection_node_end_id'].isin(valid_connection_node_ids))
        ]
        
        # Special treatment for cross section locations
        filtered_cross_section_location = cross_section_location[
            cross_section_location['channel_id'].isin(filtered_channel['id'])
        ]
        
        # Create new selection of connection nodes: only keep connection node connected to element, no 'floating' connection nodes
        valid_connection_node_ids = set()
        structures = [filtered_manhole, filtered_channel, filtered_pipe, filtered_orifice, 
                      filtered_culvert, filtered_weir, filtered_pumpstation]
        column_names = ['connection_node_id', 'connection_node_start_id', 'connection_node_end_id']
        
        for structure in structures:
            for column_name in column_names:
                if column_name in structure.columns:
                    valid_connection_node_ids.update(structure[column_name])
                    
        filtered_connection_node = connection_node[
            connection_node['id'].isin(valid_connection_node_ids)
        ]        
        
        valid_connection_node_ids = filtered_connection_node['id']
        
        # Continue for the new set of connection_nodes with the other items
        
        filtered_boundary_condition_1d = boundary_condition_1d[
            boundary_condition_1d['connection_node_id'].isin(valid_connection_node_ids)
        ]
        
        filtered_lateral_1d = lateral_1d[
            lateral_1d['connection_node_id'].isin(valid_connection_node_ids)
        ]
        
        filtered_impervious_surface_map = impervious_surface_map[
            impervious_surface_map['connection_node_id'].isin(valid_connection_node_ids)
        ]
        
        valid_impervious_surface_ids = filtered_impervious_surface_map['impervious_surface_id']
        filtered_impervious_surface = impervious_surface[
            impervious_surface['id'].isin(valid_impervious_surface_ids)
        ]    
        
        # Filter other items that should be completely within mask
        filtered_lateral_2d = self.spatial_join(lateral_2d, subarea_gdf, 'inner', 'within', '_subarea')
        filtered_boundary_condition_2d = self.spatial_join(boundary_condition_2d, subarea_gdf, 'inner', 'within', '_subarea')
        filtered_potential_breach = self.spatial_join(potential_breach, subarea_gdf, 'inner', 'within', '_subarea')
        
        # Filter other items that should be atleast partly within mask
        filtered_linear_obstacle = self.spatial_join(linear_obstacle, subarea_gdf, 'inner', 'intersects', '_subarea')
        filtered_exchange_line = self.spatial_join(exchange_line, subarea_gdf, 'inner', 'intersects', '_subarea')
        filtered_grid_refinement = self.spatial_join(grid_refinement, subarea_gdf, 'inner', 'intersects', '_subarea')
        filtered_grid_refinement_area = self.spatial_join(grid_refinement_area, subarea_gdf, 'inner', 'intersects', '_subarea')

        # Write all filtered items to the schematisation gpkg
        filtered_connection_node.to_file(output_schematisation_gpkg_path, layer='connection_node', driver="GPKG")
        filtered_manhole.to_file(output_schematisation_gpkg_path, layer='manhole', driver="GPKG")
        filtered_pipe.to_file(output_schematisation_gpkg_path, layer='pipe', driver="GPKG")
        filtered_weir.to_file(output_schematisation_gpkg_path, layer='weir', driver="GPKG")
        filtered_orifice.to_file(output_schematisation_gpkg_path, layer='orifice', driver="GPKG")
        filtered_culvert.to_file(output_schematisation_gpkg_path, layer='culvert', driver="GPKG")
        filtered_cross_section_location.to_file(output_schematisation_gpkg_path, layer='cross_section_location', driver="GPKG")
        filtered_channel.to_file(output_schematisation_gpkg_path, layer='channel', driver="GPKG")
        filtered_pumpstation_map.to_file(output_schematisation_gpkg_path, layer='pumpstation_map', driver="GPKG")
        filtered_pumpstation.to_file(output_schematisation_gpkg_path, layer='pumpstation', driver="GPKG")
        filtered_boundary_condition_1d.to_file(output_schematisation_gpkg_path, layer='1d_boundary_condition', driver="GPKG")
        filtered_boundary_condition_2d.to_file(output_schematisation_gpkg_path, layer='2d_boundary_condition', driver="GPKG")
        filtered_lateral_1d.to_file(output_schematisation_gpkg_path, layer='1d_lateral', driver="GPKG")
        filtered_lateral_2d.to_file(output_schematisation_gpkg_path, layer='2d_lateral', driver="GPKG")
        filtered_impervious_surface_map.to_file(output_schematisation_gpkg_path, layer='impervious_surface_map', driver="GPKG")
        filtered_impervious_surface.to_file(output_schematisation_gpkg_path, layer='impervious_surface', driver="GPKG")
        filtered_linear_obstacle.to_file(output_schematisation_gpkg_path, layer='linear_obstacle', driver="GPKG")
        filtered_potential_breach.to_file(output_schematisation_gpkg_path, layer='potential_breach', driver="GPKG")
        filtered_exchange_line.to_file(output_schematisation_gpkg_path, layer='exchange_line', driver="GPKG")
        filtered_grid_refinement.to_file(output_schematisation_gpkg_path, layer='grid_refinement', driver="GPKG")
        filtered_grid_refinement_area.to_file(output_schematisation_gpkg_path, layer='grid_refinement_area', driver="GPKG")
    
        ### raster clipping
        rasters = False
        for file_path in self.rasters_directory.iterdir():
            if file_path.is_file() and file_path.suffix.lower() in {'.tif', '.tiff'}:
                rasters = True
                break
                
        if rasters:          
            # Select the 'rekencellen' features that intersect with the current 'subarea'
            intersecting_calculation_grid_cells = self.calculation_grid_cells[self.calculation_grid_cells.intersects(subarea.geometry)]
        
            # Check if there are intersecting 'rekencellen'
            if intersecting_calculation_grid_cells.empty:
                raise NoCalcGridCellsSelected(f"No intersecting calculation grid cells for subarea '{name}'. Check your subarea extent.")
                
            else:
                # Union the intersecting 'rekencellen' to create a single polygon
                dissolved_geometry = unary_union(intersecting_calculation_grid_cells.geometry)
        
                for file_path in self.rasters_directory.iterdir():
                    if file_path.is_file() and file_path.suffix.lower() in {'.tif', '.tiff'}:
                        output_folder = output_directory / 'rasters' 
                        output_folder.mkdir(parents=True, exist_ok=True)
                        output_path = output_folder / file_path.name
                        self.clip_raster(file_path, dissolved_geometry, output_path)

    def read_geopackage_layers(self, gpkg_path):
        layers_dict = {}
        
        for subareas_layer_name in fiona.listlayers(gpkg_path):
            layer_gdf = gpd.read_file(gpkg_path, layer=subareas_layer_name, driver='GPKG')
            layers_dict[subareas_layer_name] = layer_gdf
        
        return layers_dict
    
    def spatial_join(self, layer, mask, how, predicate, rsuffix):
        filtered_layer = gpd.sjoin(
            layer, mask, how=how, predicate=predicate, rsuffix=rsuffix
        )
        
        clean_filtered_layer = filtered_layer[layer.columns]
        
        return clean_filtered_layer
    
    def clip_raster(self, input_path, mask_geometry, output_path):
        """Clips raster based on mask (shapely polygon geometry)."""
        
        # Open input raster and get CRS
        src_ds = gdal.Open(str(input_path))
        input_crs = osr.SpatialReference()
        input_crs.ImportFromWkt(src_ds.GetProjection())
        
        # Convert mask geometry to GeoJSON
        mask_geojson = mapping(mask_geometry)
        
        # Set the coordinate reference system for the mask geometry
        mask_geojson['crs'] = {
            'type': 'name',
            'properties': {'name': input_crs.ExportToWkt()}
        }
        
        # Create a temporary GeoJSON file for the mask geometry
        temp_geojson_file = tempfile.NamedTemporaryFile(suffix=".geojson", delete=False)
        with open(temp_geojson_file.name, 'w') as geojson_file:
            json.dump(mask_geojson, geojson_file)
        
        # Clip raster
        options = gdal.WarpOptions(
            format='GTiff',
            cutlineDSName=str(temp_geojson_file.name),
            #cutlineLayer='OGRGeoJSON',
            cropToCutline=True,
            dstSRS=input_crs.ExportToWkt(),
            creationOptions=['COMPRESS=DEFLATE']
        )
        clipped_ds = gdal.Warp(str(output_path), src_ds, options=options)
        
        # Close datasets
        src_ds = None
        clipped_ds = None
        
        # Clean up temporary GeoJSON file
        temp_geojson_file.close()

                
def run(schematisation_directory, subareas_path, field_name, calculation_grid_cells_path, 
        subareas_layer_name=None, calculation_grid_cells_layer_name=None):

    submodels = Submodels(
        schematisation_directory = schematisation_directory, 
        subareas_path = subareas_path, 
        calculation_grid_cells_path = calculation_grid_cells_path,
        field_name = field_name, 
        subareas_layer_name = subareas_layer_name, 
        calculation_grid_cells_layer_name = calculation_grid_cells_layer_name
    )
    
def get_parser():
    """Return argument parser."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "schematisation_directory", help="Required. Folder including .sqlite file, .gpkg file and rasters folder."
    )
    
    parser.add_argument(
        "subareas_path",
        help="Required. Path to vector containing geometry of submodel areas as different features. Most common shape file extensions supported.",
    )
       
    parser.add_argument(
        "field_name",
        help = "Required. Name of field containing names of submodels. This name will be used to name output folder and submodel. Should contain unique subarea names."
        )
    
    parser.add_argument(
        "calculation_grid_cells_path",
        help="Required. Path to vector containing geometry of calculation grid cells. Most common shape file extensions supported.",
    )
    
    parser.add_argument(
        "--subareas_layer_name", "-s", 
        help="Optional. Will be used as layer name in case subareas_path is of format geopackage (.gpkg)",
    )
    
    parser.add_argument(
        "--calculation_grid_cells_layer_name", "-c",
        help="Optional. Will be used as layer name in case calculation_grid_cells_path is of format geopackage (.gpkg)",
    )
    
    return parser


def main():
    """Call extract_all with args from parser."""
    return run(**vars(get_parser().parse_args()))


if __name__ == "__main__":
    exit(main())    


