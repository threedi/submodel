from pathlib import Path
from osgeo import ogr
import logging

# remove logging file if it exists
if Path("add_model_2_model.log").exists():
    Path("add_model_2_model.log").unlink()

# Configure logging to file and console
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler("add_model_2_model.log"), logging.StreamHandler()],
)


class GeoPackageFileNotFoundError(Exception):
    pass


# These layers should contain one feature only
# Assume base layer contains the valid values
BASE_ONLY_LAYERS = [
    "schema_version",
    "initial_conditions",
    "interception",
    "interflow",
    "groundwater",
    "simple_infiltration",
    "vegetation_drag_2d",
    "model_settings",
    "aggregation_settings",
    "numerical_settings",
    "physical_settings",
    "simulation_template_settings",
    "time_step_settings",
]

CREATE_LAYER_OPTIONS = ["GEOMETRY_NAME=geom", "FID=id"]

RELATIONSHIPS = {
    "channel": {
        "id": [
            ("cross_section_location", ["channel_id"]),
        ]
    },
    "connection_node": {
        "id": [
            ("boundary_condition_1d", ["connection_node_id"]),
            ("lateral_1d", ["connection_node_id"]),
            ("channel", ["connection_node_id_start", "connection_node_id_end"]),
            ("pump", ["connection_node_id"]),
            ("pump_map", ["connection_node_id_end"]),
            ("weir", ["connection_node_id_start", "connection_node_id_end"]),
            ("orifice", ["connection_node_id_start", "connection_node_id_end"]),
            ("culvert", ["connection_node_id_start", "connection_node_id_end"]),
            ("pipe", ["connection_node_id_start", "connection_node_id_end"]),
        ]
    },
}


class AddModel2Model:
    def __init__(
        self,
        base_schematisation_directory: Path,
        additional_schematisation_directory: Path,
        output_schematisation_directory: Path,
    ):
        self.base_dir = base_schematisation_directory
        self.add_dir = additional_schematisation_directory
        self.out_dir = output_schematisation_directory
        self.out_dir.mkdir(parents=True, exist_ok=True)

        self.base_gpkg = self._find_gpkg(self.base_dir)
        self.add_gpkg = self._find_gpkg(self.add_dir)

        self.base_layers = self._read_layers(self.base_gpkg)
        self.add_layers = self._read_layers(self.add_gpkg)

        self._merge_layers()
        self._save_layers()

    def _find_gpkg(self, directory: Path) -> Path:
        files = list(directory.glob("*.gpkg"))
        if not files:
            raise GeoPackageFileNotFoundError(
                f"No GeoPackage file found in {directory}."
            )
        return files[0]

    def _read_layers(self, gpkg_path: Path) -> dict:
        ds = ogr.Open(str(gpkg_path))
        layers = {}
        for i in range(ds.GetLayerCount()):
            layer = ds.GetLayerByIndex(i)
            features = [f.Clone() for f in layer]
            layers[layer.GetName()] = features
        ds = None
        return layers

    def _get_pk(self, feature):
        return feature.GetFID()

    def _merge_layers(self):
        print("Merging layers...")
        merged = {k: v.copy() for k, v in self.base_layers.items()}
        # Store all old2new mappings for each layer
        layer_old2new = {}

        for lname, add_feats in self.add_layers.items():
            if lname not in merged:
                merged[lname] = [f.Clone() for f in add_feats]
                continue

            if lname in BASE_ONLY_LAYERS:
                continue

            base_feats = merged[lname]
            max_pk = max((self._get_pk(f) for f in base_feats), default=0)
            old2new = {}
            new_feats = []
            for i, feat in enumerate(add_feats):
                old_pk = self._get_pk(feat)
                new_pk = max_pk + 1 + i
                old2new[old_pk] = new_pk
                feat.SetFID(new_pk)
                new_feats.append(feat.Clone())

            if new_feats:
                merged[lname] = base_feats + new_feats
            else:
                merged[lname] = base_feats

            # Save mapping for this layer
            layer_old2new[lname] = old2new

            logging.info(
                f"Merged layer '{lname}': {len(base_feats)} base features, {len(add_feats)} added features, {len(merged[lname])} total features."
            )

        # After all layers are merged update relationships based on old2new mappings
        # For instance, we need to update 'channel_id' in 'cross_section_location'.
        for lname, old2new in layer_old2new.items():
            if lname not in RELATIONSHIPS.keys():
                continue

            logging.info(f"Processing relationships for layer '{lname}'.")

            for parent_field, mapper in RELATIONSHIPS[lname].items():
                for rel_layer, rel_fields in mapper:
                    if rel_layer not in merged:
                        logging.warning(
                            f"Relationship layer '{rel_layer}' not found in merged layers."
                        )
                        continue

                    rel_feats = merged[rel_layer]

                    logging.info(
                        f"Updating relationships in layer '{rel_layer}' for fields {rel_fields} based on '{lname}.{parent_field}'"
                    )

                    count_updated = 0
                    count_not_found = 0
                    for rel_feat in rel_feats:
                        for rel_field in rel_fields:
                            if rel_feat.IsFieldSet(rel_field):
                                old_value = rel_feat.GetField(rel_field)
                                if old_value in old2new:
                                    new_value = old2new[old_value]
                                    rel_feat.SetField(rel_field, new_value)
                                    count_updated += 1
                                else:
                                    count_not_found += 1

                    logging.info(
                        f"Updated {count_updated} features in layer '{rel_layer}' for fields {rel_fields}. {count_not_found} values not found in old2new mapping."
                    )

        self.merged_layers = merged

    def _save_layers(self):
        print("Saving merged layers to GeoPackage...")
        out_gpkg = self.out_dir / "merged_schematisation.gpkg"
        driver = ogr.GetDriverByName("GPKG")
        if out_gpkg.exists():
            driver.DeleteDataSource(str(out_gpkg))
        ds = driver.CreateDataSource(str(out_gpkg))

        for lname, feats in self.merged_layers.items():
            if not feats:
                print(f"No features to write for layer {lname}. Copying from base.")
                copy_layer_between_geopackages(
                    str(self.base_gpkg), lname, str(out_gpkg), lname
                )
                continue

            # Determine if the layer is spatial or non-spatial
            first_geom_feat = next(
                (f for f in feats if f.GetGeometryRef() is not None), None
            )
            if first_geom_feat is not None:
                geom_type = first_geom_feat.GetGeometryRef().GetGeometryType()
                srs = (
                    first_geom_feat.GetGeometryRef().GetSpatialReference()
                    if first_geom_feat.GetGeometryRef()
                    else None
                )
            else:
                # Non-spatial layer
                geom_type = ogr.wkbNone
                srs = None

            base_layer_feats = self.base_layers.get(lname, [])
            if not base_layer_feats:
                continue

            base_defn = base_layer_feats[0].GetDefnRef()
            field_defs = {
                base_defn.GetFieldDefn(i).GetName(): base_defn.GetFieldDefn(i)
                for i in range(base_defn.GetFieldCount())
            }

            layer = ds.CreateLayer(
                lname, srs=srs, geom_type=geom_type, options=CREATE_LAYER_OPTIONS
            )

            for fname, fdef in field_defs.items():
                fd = ogr.FieldDefn(fname, fdef.GetType())
                layer.CreateField(fd)

            # Start transaction for this layer
            layer.StartTransaction()
            try:
                layer_defn = layer.GetLayerDefn()
                for feat in feats:
                    out_feat = ogr.Feature(layer_defn)
                    for fname in field_defs:
                        if feat.IsFieldSet(fname):
                            out_feat.SetField(fname, feat.GetField(fname))
                    out_feat.SetFID(feat.GetFID())
                    geom = feat.GetGeometryRef()
                    if geom is not None and geom_type != ogr.wkbNone:
                        out_feat.SetGeometry(geom.Clone())
                    layer.CreateFeature(out_feat)
                    out_feat = None
                layer.CommitTransaction()
            except Exception as e:
                layer.RollbackTransaction()
                print(f"Error writing layer {lname}: {e}")
        ds = None


def copy_layer_between_geopackages(
    source_gpkg: str, source_layer_name: str, target_gpkg: str, target_layer_name: str
) -> None:
    # Open source GPKG (read-only)
    src_ds = ogr.Open(source_gpkg, 0)
    if not src_ds:
        raise RuntimeError(f"Failed to open source: {source_gpkg}")

    src_layer = src_ds.GetLayerByName(source_layer_name)
    if not src_layer:
        raise ValueError(f"Layer '{source_layer_name}' not found in {source_gpkg}")

    # Open target GPKG (read/write)
    tgt_ds = ogr.Open(target_gpkg, 1)
    if not tgt_ds:
        raise RuntimeError(f"Failed to open target: {target_gpkg}")

    # Check if target layer already exists
    if tgt_ds.GetLayerByName(target_layer_name):
        raise ValueError(f"Layer '{target_layer_name}' already exists in {target_gpkg}")

    # Create new layer
    tgt_layer = tgt_ds.CreateLayer(
        target_layer_name,
        geom_type=src_layer.GetGeomType(),
        srs=src_layer.GetSpatialRef(),
        options=CREATE_LAYER_OPTIONS,
    )

    # Copy fields
    src_defn = src_layer.GetLayerDefn()
    for i in range(src_defn.GetFieldCount()):
        field_defn = src_defn.GetFieldDefn(i)
        tgt_layer.CreateField(field_defn)

    # Copy features
    for feature in src_layer:
        new_feature = ogr.Feature(tgt_layer.GetLayerDefn())
        new_feature.SetFrom(feature)
        tgt_layer.CreateFeature(new_feature)
        new_feature = None

    # Clean up
    src_ds = None
    tgt_ds = None


# Example usage:
if __name__ == "__main__":
    base_dir = Path(
        r"C:\Users\stijn.overmeen\testing\Y088 Sallandse Wetering Hoofdmodel\work in progress\schematisation"
    )
    add_dir = Path(
        r"C:\Users\stijn.overmeen\testing\Y088 Sallandse Wetering BC Grof Inputmodel\work in progress\schematisation\gebied_noordkant"
    )
    out_dir = Path(
        r"C:\Users\stijn.overmeen\testing\Y088 Sallandse Wetering BC Grof Inputmodel\work in progress\schematisation\merge"
    )

    try:
        merger = AddModel2Model(base_dir, add_dir, out_dir)
        print("Schematisations merged successfully.")
    except GeoPackageFileNotFoundError as e:
        print(e)
