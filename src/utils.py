import ee
import geemap
import os
import pandas as pd
from datetime import datetime
from google.cloud import storage


class GeometryOperations:
    def __init__(self, buffer_distance=10000, max_error=1):
        self.buffer_distance = buffer_distance
        self.max_error = max_error
        self.water_mask = ee.Image('JRC/GSW1_0/GlobalSurfaceWater')

    def buffer_polygon(self, feat):
        """Create buffer around polygon"""
        feat = ee.Feature(feat)
        out = feat.buffer(self.buffer_distance).geometry()
        inn = feat.buffer(-self.buffer_distance).geometry()
        aoi = out.difference(inn, self.max_error)
        return aoi

    def mask_water(self, feat):
        """Mask water bodies from feature"""
        water_no_holes = self.water_mask.select('max_extent')\
            .focalMax(radius=30, units='meters', kernelType='square')\
            .focalMin(radius=30, units='meters', kernelType='square')
        water_vect = water_no_holes.reduceToVectors(
            reducer=ee.Reducer.countEvery(),
            geometry=feat.buffer(1000),
            scale=30,
            maxPixels=1e10,
            geometryType='polygon',
            eightConnected=False)
        geom = feat.difference(water_vect.geometry(), maxError=self.max_error)
        return geom
    
    def get_biome(self, geom): 
        """Get biome with largest overlap for a feature, add BIOME_NAME property"""
        ecoregions = ee.FeatureCollection("RESOLVE/ECOREGIONS/2017")
        intersecting = ecoregions.map(lambda eco: eco.set(
            'intersection_area', 
            eco.geometry().intersection(geom).area()
        )).filterBounds(geom)
        largest_ecoregion = intersecting.sort('intersection_area', False).first()
        result = ee.Feature(geom).set('BIOME_NAME', largest_ecoregion.get('BIOME_NAME'))
        return result
    
    def get_pixels_boundary(self, image, polygon, scale=10):
        """Get pixels that straddle the 10km surrounding the polygon boundary"""
        outer_buffer = ee.Geometry(polygon).buffer(scale/2)
        inner_buffer = ee.Geometry(polygon).buffer(-scale/2)
        boundary_region = outer_buffer.difference(inner_buffer)
        boundary_pixels = image.updateMask(
            image.clip(boundary_region)
            .mask()
            .reduce(ee.Reducer.anyNonZero())
        )
        return boundary_pixels
        

class ImageOperations:
    def __init__(self):
        self.modis = ee.ImageCollection('MODIS/006/MOD09A1')

    def filter_for_year(self, feat, year):
        """Filter images for specific year"""
        start = ee.Date.fromYMD(year, 1, 1)
        return ee.Filter.And(
            ee.Filter.bounds(feat),
            ee.Filter.date(start, start.advance(1, "year"))
        )

    def add_indices_to_image(self, image):
        """Add vegetation indices to image"""
        EVI = image.expression(
            "2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))",
            {
                'NIR': image.select("sur_refl_b02"),
                'RED': image.select("sur_refl_b01"),
                'BLUE': image.select("sur_refl_b03")
            }
        ).rename("EVI")

        NDVI = image.expression(
            "(NIR - RED) / (NIR + RED)",
            {
                'NIR': image.select("sur_refl_b02"),
                'RED': image.select("sur_refl_b01")
            }
        ).rename("NDVI")

        return image.addBands([EVI, NDVI])

    def get_gradient_magnitude(self, image):
        """Calculate gradient magnitude"""
        gradient = image.gradient()
        gradient_x = gradient.select('x')
        gradient_y = gradient.select('y')
        magnitude = gradient_x.pow(2).add(gradient_y.pow(2)).sqrt()
        return magnitude


class StatsOperations:
    def __init__(self):
        self.gHM_collection = ee.ImageCollection('CSP/HM/GlobalHumanModification')

    def calculate_gradient_statistics(self, layer, scale=500, name='buffer'):
        """Calculate mean and standard deviation of gradient magnitude"""
        bounded_layer = layer.clip(layer.geometry())
        stats = bounded_layer.reduceRegion(
            reducer=ee.Reducer.mean().combine(
                reducer2=ee.Reducer.stdDev(),
                sharedInputs=True
            ).combine(
                reducer2=ee.Reducer.sum(),
                sharedInputs=True
            ),
            geometry=bounded_layer.geometry(),
            scale=scale,
            maxPixels=1e10
        )
        return stats

    def calculate_human_modification(self, geometry, scale=500):
        """Calculate mean Global Human Modification value"""
        mean_gHM = self.gHM_collection.mean().clip(geometry)
        gHM_value = mean_gHM.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=geometry,
            scale=scale,
            maxPixels=1e9)
        return gHM_value.get('gHM')


class FeatureProcessor:
    def __init__(self, geo_ops, img_ops, stats_ops):
        self.geo_ops = geo_ops
        self.img_ops = img_ops
        self.stats_ops = stats_ops
        self.bands_to_process = ['sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03', 'EVI', 'NDVI']
        
    def collect_feature_info(self, pa, aoi_with_biome):
        """Collect basic protected area feature information"""
        return {
            'WDPA_PID': pa.get('WDPA_PID'),
            'ORIG_NAME': pa.get('ORIG_NAME'),
            'BIOME_NAME': aoi_with_biome.get('BIOME_NAME'),
            'GIS_AREA': pa.get('GIS_AREA')
        }
    
    def process_single_band(self, band_name, image, pa_geometry):
        """Process a single band and return its statistics"""
        single_band = image.select(band_name)
        gradient = self.img_ops.get_gradient_magnitude(single_band)
        boundary_pixels = self.geo_ops.get_pixels_boundary(gradient, pa_geometry, scale=500)
        boundary_pixels = boundary_pixels.clip(pa_geometry.buffer(500))
        return {
            'band_name': band_name,
            'boundary_stats': self.stats_ops.calculate_gradient_statistics(boundary_pixels, name='boundary'),
            'buffer_stats': self.stats_ops.calculate_gradient_statistics(boundary_pixels, name='buffer'),
            'gradient': gradient,
            'boundary_pixels': boundary_pixels}
    
    def process_all_bands(self, image, pa_geometry):
        """Process all bands and collect their statistics"""
        return [self.process_single_band(band_name, image, pa_geometry) 
                for band_name in self.bands_to_process]
    
    def compile_statistics(self, feature_info, computed_stats, year):
        """Compile all statistics into a list of dictionaries"""
        feature_values = {k: v.getInfo() for k, v in feature_info.items()}
        all_stats = []
        for stat in computed_stats:
            row_stats = {
                **feature_values,
                'band_name': stat['band_name'],
                'year': year,
                **{f"boundary_{k}": v for k, v in stat['boundary_stats'].getInfo().items()},
                **{f"buffer_{k}": v for k, v in stat['buffer_stats'].getInfo().items()}}
            all_stats.append(row_stats)       
        return all_stats


class ExportResults: 
    def __init__(self, results_path='/workspaces/jupyter-gee-project/results'):
        self.results_path = results_path
        os.makedirs(results_path, exist_ok=True)
    
    def generate_filename(self, protected_area_name, year):
        """Generate standardized filename for results."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return f"results_{protected_area_name.replace(' ', '_')}_{year}_{timestamp}.csv"
    
    def save_statistics_to_csv(self, all_stats, protected_area_name, year):
        """Save computed statistics to a CSV file"""
        df = pd.DataFrame(all_stats)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_file = f'results_{protected_area_name.replace(" ", "_")}_{year}_{timestamp}.csv'
        df.to_csv(output_file, index=False)
        print(f'Results saved to {output_file}')
        return df, output_file

    def save_df_to_gcs(self, df, bucket_name, wdpaid, year):
        """Save DataFrame as CSV and upload to GCS."""
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        tmp_file = "/workspace/temp.csv"
        df.to_csv(tmp_file, index=False)
        client = storage.Client(project=bucket_name)
        bucket = client.bucket(bucket_name)
        blob = bucket.blob(f'protected_areas/tables/_{wdpaid.replace(" ", "_")}_{year}_{timestamp}.csv')
        blob.upload_from_filename(tmp_file)
        print(f"Uploaded to: gs://{bucket_name}/{blob.name}")

    
    def export_image_to_cloud(image: ee.Image, date: str):
        path_file = 'landsat_lebanon/landsat_indices_' + date
        print(path_file)
        export_task = ee.batch.Export.image.toCloudStorage(
        image=image,
        description='landsat_export_cog',
        bucket='dse-staff', 
        fileNamePrefix=path_file,  
        fileFormat='GeoTIFF', 
        formatOptions={
            'cloudOptimized': True,  
        },
        maxPixels=1e8,  
        scale=30  
        )
        export_task.start()
        return print(f'{date} saved')


class Visualization:
    def __init__(self):
        self.default_vis_params = {
            'min': -0.5,
            'max': 1,
            'palette': ['black', 'gray', 'white']
        }

    def create_map(self, geometry, gradient, boundary_pixels, vis_params=None):
        """Create and return an interactive map"""
        if vis_params is None:
            vis_params = self.default_vis_params
        Map = geemap.Map()
        Map.add_basemap('HYBRID')
        Map.centerObject(geometry, 8)
        Map.addLayer(geometry, {'color': 'red'}, 'Protected Area Geometry')
        Map.addLayer(gradient, vis_params, 'Gradient Layer')
        Map.addLayer(boundary_pixels, vis_params, 'Gradient Boundary Pixels')
        return Map