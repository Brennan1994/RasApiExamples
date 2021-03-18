using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using H5Assist;
using Geospatial.Rasters;
using Geospatial.GDALAssist;
using RasMapperLib;
using RasMapperLib.Mapping;
using RasMapperLib.Utilities;
using Utility.Progress;
using Path = System.IO.Path;

/* This project requires adding references to System.Windows.Forms because of RasMapperLib. Apologies.
 * Easiest way to get it going is to install a RAS v6 Beta, which will install the necessary pre-reqs. If you want to scrape
 * the files and install prereqs manually, they are...
 * .Net 4.7.2 Runtime   (https://dotnet.microsoft.com/download/dotnet-framework/net472)
 * VC Runtime 2015-2019 (https://support.microsoft.com/en-us/topic/the-latest-supported-visual-c-downloads-2647da03-1eea-4433-9aff-95f26a218cc0) (x86/x64 as needed)
 * 
 * Shipped with necessary references and a gdal build.
 * 
 * Anywhere you see "RasMapperLib.Utilities.Converter.Convert(blah)", that's transitioning between new/old world equivalent code.
 * 
 * RasMapperLib is old world, Geospatial is new world.
 */
namespace RasApiExamples
{
    class Program
    {
        const string OutputDirectory = @"C:\Temp\RasOutput\";

        static void Main(string[] args)
        {
            // Point this to your RAS install directory
            string installDirectory = ".\\";
            string gdalDirectory = Path.Combine(installDirectory, "GDAL");

            if (!Directory.Exists(gdalDirectory))
            {
                Console.WriteLine("GDAL directory not found: " + gdalDirectory);
                return;
            }

            // Initializing GDAL is required for many common RasMapper operations. 
            // The InitializeMultiplatform call is new (and required) for any of the Ras version 6 builds.
            // To target an older v5 distribution, you'll have to use the "Initialize" call, and point it
            // manually to the x86/x64 sub-directory.
            GDALSetup.InitializeMultiplatform(gdalDirectory);

            string resDirectory = @"C:\Users\q0hecbbb\Projects\FDA vs Event Tree LifeSim\Data Collection\West_Sac_Material\Hydraulics\RAS\West_Sac_HECRAS\RAS_Model";
            var resFiles = Directory.GetFiles(resDirectory, "*.p*.hdf");

            string terrFilez = @"C:\Users\q0hecbbb\Projects\FDA vs Event Tree LifeSim\Data Collection\West_Sac_Material\Hydraulics\RAS\West_Sac_HECRAS\RAS_Model\Terrain\Terrain.hdf";
            string ptShpFile = @"C:\Users\q0hecbbb\Projects\FDA vs Event Tree LifeSim\Working\GIS\Final_Structure_Inventory.shp";
            string projShpPrj = @"C:\Users\q0hecbbb\Projects\FDA vs Event Tree LifeSim\Working\RAS\West_Sac_HECRAS\RAS_Model - Copy\Terrain\NAD 1983 StatePlane California II FIPS 0402 (US Feet).prj";

            SampleGridExportShapefile.ExportMaxDepthShapefile(resFiles, terrFilez, ptShpFile,projShpPrj, @"C:\Temp\ForBrennanOutput.shp");
            return;

            // Sample result/terrain.
            string resultFile = @"C:\Work\_datasets\_core\Muncie\Muncie.p03.hdf";
            string terrFile = @"C:\Work\_datasets\_core\Muncie\Terrain\TerrainWithChannel.hdf";
            string projectionFile = @"C:\Work\_datasets\_core\Muncie\Terrain\NAD 1983 StatePlane Indiana East FIPS 1301 (US Feet).prj";

            // Sample Project
            // string rasmapFile = @"C:\Work\_datasets\_core\Muncie\Muncie.rasmap";
            // SummarizeProject(rasmapFile);


            // Generate a stored depth map from the default terrain association
            var res = new RASResults(resultFile);
            GenerateStoredMap(res, Path.Combine(OutputDirectory, "MaxDepth"), MapTypes.Depth, RASResultsMap.MaxProfileIndex);


            // Generate a stored depth map with a different terrain. (Also necessary if the terrain has been moved)

            // Note the first argument is the terrain layer name, which doesn't mean much in this context.
            // Constructor signature likely to change for v7.
            var terrain = new TerrainLayer("Terrain", terrFile);
            GenerateStoredMap(res, terrain, Path.Combine(OutputDirectory, "MaxDepth_SpecifiedTerrain"), MapTypes.Depth, RASResultsMap.MaxProfileIndex);



            // Dynamic map query (raster)
            var map = new RASResultsMap(res, MapTypes.Depth);
            map.ProfileIndex = RASResultsMap.MaxProfileIndex;
            var geom = res.Geometry;
            var terr = geom.Terrain;

            // Get the extents of the geometry
            var ext = geom.Extent;

            // Create a raster/grid definition with a coarse resolution
            // 10k pixels is ~100x100 if the domain is square.
            var rasterDef = RasterM.ComputeRasterM(10000, ext);

            float[] rasterData = null;
            map.ResampleDynamic(rasterDef, ref rasterData);

            // Note: We always use InterpolatedLayer.NoData for internal no-data values. (-9999)
            // Any dry pixels will be set to this value.


            // Dynamic map query (points)
            var center = ext.CenterPointM();
            var ptList = new PointMs();
            ptList.Add(center);
            ptList.Add(center.Translate(-100, 0));
            ptList.Add(center.Translate(100, 0));
            ptList.Add(center.Translate(0, -100));
            ptList.Add(center.Translate(0, 100));

            // Sample the given points to see where they intersect the geometry (XS/SA/2D)
            var mapPixels = geom.MapPixels(ptList);

            // Sample the terrain at this point-list
            var terrainElevations = terr.ComputePointElevations(ptList);

            // Compute the map values at these point locations. Build an HDF5 CacheCollection
            // if you're doing many consecutive calls on the same file. (It's set to null here)
            float[] pointData = null;
            res.ComputeSwitch(map, mapPixels, RASResultsMap.MaxProfileIndex, terrainElevations, null, ref pointData);



            // Query time-series data. This will create a secondary file called PostProcessing.hdf that is path-relative
            // beneath the result file, e.g. ./[Result Name]/PostProcessing.hdf.
            // This file does a vertical striation of the result data, making time-series calls dramatically faster.
            // If this isn't acceptable, you can always loop through every timestep and invoke the above compute call.
            // That can be several orders of magnitude slower for large datasets.
            float[][] pointTimeSeries = null;

            // Returned buffer is float[PointCount][ProfileCount]
            res.ComputeSwitchTimeSeries(map, mapPixels, terrainElevations, null, ref pointTimeSeries, false, false);


            // Export to an image
            var sfill = MapTypes.Depth.DefaultSurfaceFill;

            // If you just want colored pixels (ARGB)
            //int[] colorBuf = new int[rasterData.Length];
            //sfill.ParallelPaint(InterpolatedLayer.NoData, rasterData, colorBuf);

            string imgFile = Path.Combine(OutputDirectory, "image.png");
            using (var bmp = sfill.ComputeBitmap(InterpolatedLayer.NoData, rasterData, rasterDef.Cols, rasterDef.Rows))
            {
                bmp.Save(imgFile, System.Drawing.Imaging.ImageFormat.Png);
            }


            // Write our point-value depths out to a shapefile
            string shpFile = Path.Combine(OutputDirectory, "points.shp");
            var ptlyr = new PointFeatureLayer();
            for (int i = 0; i < ptList.Count; i++)
            {
                var ptm = ptList[i];
                ptm.Z = pointData[i]; // Assign the max depth

                ptlyr.AddFeature(new Point(ptm));
            }

            ptlyr.SourceFilename = shpFile;
            ptlyr.Save();


            // Read and print
            var readPointLayer = new PointFeatureLayer("points", shpFile);
            for (int i = 0; i < readPointLayer.FeatureCount(); i++)
            {
                var pt = readPointLayer.Point(i);
                Console.WriteLine("Point " + i.ToString() + ": " + pt.ToString());
            }


            SampleReprojectionToWGS84(projectionFile, map);
        }

        static void SampleReprojectionToWGS84(string srcProjectionFile, RASResultsMap map)
        {
            // Warp output data into WGS84
            var srcProjection = new ESRIProjection(srcProjectionFile);
            var dstProjection = new EPSGProjection(4326); // WGS84

            // Native geometry extent
            Extent ext = map.Results.Geometry.Extent;

            // Warp native extent to wgs84. This samples multiple points along the perimeter of the given extent and returns the bounding box.
            var wgs84Ext = Converter.Convert(Projector.GetBoundingBox(Converter.Convert(ext), srcProjection, dstProjection));

            // Can target a specific cell-size, or number of cells.
            const double ArcSecond = 1d / 60d / 60d; // 1 arc-second, ~100ft

            // ~10 ft cell-size
            double wgs84CellSize = ArcSecond / 10;
            int cols = (int)(wgs84Ext.Width / wgs84CellSize);
            int rows = (int)(wgs84Ext.Height / wgs84CellSize);

            const int TileSize = 256;
            int tilesWide = cols / TileSize;
            if (cols % TileSize > 0)
                tilesWide += 1;

            int tilesTall = rows / TileSize;
            if (rows % TileSize > 0)
                tilesTall += 1;

            // Raster with 1 cell per 256x256 tile
            var fullTileRaster = new RasterM(wgs84Ext.UpperLeft, tilesTall, tilesWide, wgs84CellSize * 256);

            // Raster with 1 cell per pixel (will have different extents if not perfect multiple of 256)
            // var fullPixelRaster = new RasterM(wgs84Ext.UpperLeft, rows, cols, wgs84CellSize);

            float[] nativeBuffer = null;
            float[] wgs84Buffer = new float[TileSize * TileSize];

            // Loop over each tile in wgs84 space. Warp each tile back to native projection, sample the result,
            // and then warp the data to wgs84. This could be improved by looping over larger groups of data
            // (e.g. 2x2 tiles at a time)
            for (int i = 0; i < fullTileRaster.CellCount; i++)
            {
                var tileExtent = fullTileRaster.CellExtent(i);
                var tileRaster = new RasterM(TileSize, TileSize, tileExtent);

                // Project the wgs84 tile extent back to our native projection so we can sample the map
                var srcTileExtent = Converter.Convert(Projector.GetBoundingBox(Converter.Convert(tileExtent), dstProjection, srcProjection));

                // Going across nonlinear reprojections will warp the data heavily. Getting the bounding box above (in native projection)
                // covers more area than the wgs84 tile, since the edges are effectively "curved", and we're taking the bounding box of those 
                // curved edges. To do well with the 256x256 wgs84 resolution, we should slightly oversample the native resolution. This is
                // obviously pretty hand-wavy. Might be worth some future metric of "how non-linear is this transformation", and some oversampling
                // ratio based on that.
                const int OversampleTileSize = 300;
                var srcTileRaster = new RasterM(OversampleTileSize, OversampleTileSize, srcTileExtent);

                map.ResampleDynamic(srcTileRaster, ref nativeBuffer);

                // Warp from native back to wgs84
                InterpolatedLayer.WarpData(srcProjection, srcTileRaster, nativeBuffer, dstProjection, tileRaster, wgs84Buffer);

                // Do something with the wgs84 buffer of data
                var stats = Utility.Stats.Compute(wgs84Buffer, InterpolatedLayer.NoData);
                Console.WriteLine("Tile " + (i + 1).ToString() + " / " + fullTileRaster.CellCount.ToString() +
                  " - Min=" + stats.Min.ToString() + ", Max=" + stats.Max.ToString() + ", FracNoData=" + stats.FractionNoData.ToString());
            }

        }


        static void SummarizeProject(string rasmapFile)
        {
            LoadRasMapFile(rasmapFile);
            PrintSummaryInfo();

            // No extension, will be combined with your terrain filenames
            string outputDepthMap = Path.Combine(OutputDirectory, "MaxDepth");
            string outputDB = Path.Combine(OutputDirectory, "outputDB.db");


            var sampleResult = SharedData.RasMapper.ResultsGroup.FindAllResults().FirstOrDefault(res => res.RanSuccessfully);
            if (sampleResult == null)
            {
                Console.WriteLine("Cannot generate map - no results have been successfully ran.");
            }
            else
            {
                // Alternatively, you can loop through sampleResult.ProfileNames and find the index you want to map.
                // Works with any of the given MapTypes, although some (like arrival time) don't need an index.
                GenerateStoredMap(sampleResult, outputDepthMap, MapTypes.Depth, RASResultsMap.MaxProfileIndex);

                // 12 is a moderately zoomed-out view. Compute the DB with the first 10% of profile indexes
                GenerateMapTiles(sampleResult, outputDB, MapTypes.Depth, Enumerable.Range(0, sampleResult.ProfileCount / 10).ToList(), 12);
            }
        }

        static void LoadRasMapFile(string rasmapFile)
        {
            // Set the rasmapper context.

            // This sets the static SharedData.RasMapper instance, effectively opening the project
            // so you can explore it. This is unnecessary if you know the exact files you need to work with.
            // We've tried to remove most internal references to the SharedData.RasMapper singleton, but
            // that work isn't complete yet. We've carved out safe paths for the most common operations, but
            // exploring other functionality in the API may run into issues if this isn't set.
            SharedData.RasMapFilename = rasmapFile;
            SharedData.RasMapper = new RASMapper();
            SharedData.RasMapper.LoadRASMapFile();
        }

        /// <summary>
        /// This prints summary info for associated ras project instance. Must call LoadRasMapFile before this
        /// to set the context.
        /// </summary>
        static void PrintSummaryInfo()
        {
            Console.WriteLine(" --- GEOMETRIES ---");
            var geoms = SharedData.RasMapper.GeometriesGroup.FindAllGeometries();
            foreach (var geom in geoms)
            {
                Console.WriteLine("Geometry: " + geom.Name + " (" + Path.GetFileName(geom.SourceFilename) + ")");
            }

            Console.WriteLine(Environment.NewLine + " --- RESULTS ---");
            var results = SharedData.RasMapper.ResultsGroup.FindAllResults();
            foreach (var res in results)
            {
                Console.WriteLine("Result: " + res.Name + " (" + Path.GetFileName(res.SourceFilename) + ")");

                if (!res.RanSuccessfully)
                {
                    Console.WriteLine("  (Result has not been run.)");
                    continue;
                }
                if (res.ProfileNames != null)
                {
                    Console.WriteLine("  Profiles: " + res.ProfileNames.First() + " -> " + res.ProfileNames.Last() +
                      " (" + res.ProfileCount.ToString() + " total)");
                }
            }


            Console.WriteLine(Environment.NewLine + " --- TERRAINS ---");

            var terrains = SharedData.RasMapper.TerrainsGroup.FindAllTerrains();
            foreach (var terr in terrains)
            {
                Console.WriteLine("Terrain: " + terr.Name + " (" + Path.GetFileName(terr.SourceFilename) + ")");
            }

        }

        /// <summary>
        /// Generate a stored map from the given result. The terrain must have been associated already via RasMapper,
        /// which stores the path-relative location. If the terrain has been moved since it was originally associated,
        /// you'll have to specify it manually.
        /// </summary>
        /// <param name="outputFilename">
        /// Output filename *base*, e.g. "C:\Temp\MyOutput". Any extension will be stripped. This is because we need 
        /// to create (potentially) multiple output files, one for each matching terrain file. The output files will
        /// be "C:\Temp\MyOutput.TerrainFile1.tif", "C:\Temp\MyOutput.TerrainFile2.tif", etc.
        /// </param>
        /// <param name="profileIndex">
        /// Mapping profile index to use. Can use 0 through (<see cref="RASResults.ProfileCount"/> - 1), or <see cref="RASResultsMap.MaxProfileIndex"/>
        /// or <see cref="RASResultsMap.MinProfileIndex"/>
        /// </param>
        static void GenerateStoredMap(RASResults result, string outputFilename, MapTypes mapType, int profileIndex)
        {
            var map = new RASResultsMap(result, mapType);

            // We want a stored map, with the currently associated terrain
            map.OutputMode = OutputModes.StoredDefaultTerrain;
            map.OverwriteOutputFilename = outputFilename;
            map.ProfileIndex = profileIndex;

            map.StoreMap(new ConsoleProgressReporter(true), false);
        }

        /// <summary>
        /// Generate a stored map from the given result. Forces it to use the given terrain layer instead of the
        /// internal association. Useful if the terrain has been moved, or if you're making maps against lower
        /// resolution terrains.
        /// </summary>
        static void GenerateStoredMap(RASResults result, TerrainLayer terrain, string outputFilename, MapTypes mapType, int profileIndex)
        {
            var map = new RASResultsMap(result, mapType);

            // Need to tell the map to use a different terrain than the result association
            map.OutputMode = OutputModes.StoredSpecifiedTerrain;

            // Only available in 6.0 Beta 3 or later. In prior versions, you need to set map.OverwriteTerrainLayerName, and 
            // you need RasMapper to be instantiated with the project "loaded"
            map.Terrain = terrain;

            map.OverwriteOutputFilename = outputFilename;
            map.ProfileIndex = profileIndex;

            map.StoreMap(new ConsoleProgressReporter(true), false);
        }


        /// <summary>
        /// Generates our SQLITE DB of tiles from the given result.
        /// </summary>
        static void GenerateMapTiles(RASResults result, string outputFilename, MapTypes mapType, List<int> profileIndexes, int maxZoom)
        {
            var map = new RASResultsMap(result, mapType);
            var options = new TileCacheComputeOptions(map, outputFilename, maxZoom, mapType.DefaultLayerName, profileIndexes);
            TileCacheComputable.Compute(map, options, ProgressReporter.ConsoleWrite());
        }

    }

}