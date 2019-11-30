// Google Earth Engine (GEE) script to obtain average gridMET temperatures across the HJA.
// To run it, you'll need to obtain a GEE account and upload the HJA boundaries shapefile
// in lat/long format.

var GRIDMET = ee.ImageCollection("IDAHO_EPSCOR/GRIDMET").select(['tmmn','tmmx']);// .

var Andrews = ee.FeatureCollection("users/ChrisGrayWolf/Andrews/boundaries_WGS84");

var results = GRIDMET.map(function(image) {
  return image.reduceRegions({reducer:ee.Reducer.mean(),collection:Andrews,scale:100}).map(function(f) {
    return f.setGeometry().set('date', image.date().format("YYYY-MM-dd"));
  });
});

Export.table.toDrive({collection:results.flatten(), description:'GRIDMET_reduce', folder:'Andrews'}); // ~25 min
