using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Geospatial.GDALAssist;

using RasMapperLib;
using RasMapperLib.Mapping;

namespace RasApiExamples
{
    public static class SampleGridExportShapefile
    {
        public static void ExportMaxDepthShapefile(IList<string> resultFiles, string terrainFile, string inputPrj, string inputPointShp, string outputShp)
        {
            // TODO - work with reprojection here.
            var ptLayer = new PointFeatureLayer("name", inputPointShp);

            var pts2 = ptLayer.Points().Select(p => p.PointM()).ToArray();
            var sourceProj = ShapefileStorage.GetProjection(inputPointShp);
            if(sourceProj != null)
            {
                var desproj = new ESRIProjection(inputPrj);
                var convertedPts = RasMapperLib.Utilities.Converter.Convert(pts2); //Transfers to new world to work with 7.0 Gdal Raster stuff
                Projector.TransformPoints(sourceProj, desproj, convertedPts);
                for (int i =0; i<pts2.Length; i++)
                {
                    pts2[i] = RasMapperLib.Utilities.Converter.ConvertPtM(convertedPts[i]);
                }
            }


            var terr = new TerrainLayer("terr", terrainFile);

            PointMs pts = new PointMs(pts2);
            float[] ptElevs = terr.ComputePointElevations(pts);

            // Set up the output shapefile

            ptLayer.AddAttributeColumn("Terrain", typeof(float));

            // Add every terrain elevation to the correct row in the output shapefile table
            for (int i = 0; i < ptElevs.Length; i++)
            {
                ptLayer.FeatureRow(i)["Terrain"] = ptElevs[i];
            }

            foreach (var resFn in resultFiles)
            {
                // Construct a result from the given filename.
                var res = new RASResults(resFn);
                var geo = res.Geometry;
                var wsMap = new RASResultsMap(res, MapTypes.Elevation);

                // Sample the geometry for the given points loaded from the shapefile.
                // If the geometry is the same for all of the results, we can actually reuse this object.
                // (It's pretty fast to recompute though, so I wouldn't bother)
                RASGeometryMapPoints mapPixels = geo.MapPixels(pts);

                // This will produce -9999 for NoData values.
                float[] wsValues = null;
                res.ComputeSwitch(wsMap, mapPixels, RASResultsMap.MaxProfileIndex, ptElevs, null, ref wsValues);

                // TODO - check column names.
                // I think our shapefile writer truncates to SHP character count, but this may not be unique.
                string colName = res.Name;
                ptLayer.AddAttributeColumn(colName, typeof(float));

                // Fill in the column of elevations
                for (int i = 0; i < ptElevs.Length; i++)
                {
                    ptLayer.FeatureRow(i)[colName] = wsValues[i];
                }
            }
            ptLayer.SourceFilename = outputShp;
            ptLayer.Save();
        }
    }
}
