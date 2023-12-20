 /********************************************************************************************
Take Baekdu Mountain Zone (BM) in 2021 as an example
*********************************************************************************************/ 

// Define image and table variables
var image = ee.Image("landuse_2020"), // 30m*30m land cover data (http://www.resdc.cn/)
    table = ee.FeatureCollection("China_dongbei"), // Administrative boundaries of agricultural areas in Northeast China
    table2 = ee.FeatureCollection("BM_labels_2021"); // Labels have been produced for crops
    
// Set the region of interest   
var BM_Point =ee.Geometry.Point([127.5477075829258,42.91374280135512]);
var geometry = table.filterBounds(BM_Point);


 /********************************************************************************************
Function
*********************************************************************************************/ 
// Function to rename Sentinel-2 bands, set properties
function sentinel2sr (image){
  return image.select(
                      ['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9', 'B11','B12','QA60','SCL']).rename
                    (['aerosol','blue','green','red','red1','red2','red3','nir','red4','cirrus','swir1','swir2','QA60','SCL']
                    )
                    .divide(10000).toDouble()
                    .set('solar_azimuth',image.get('MEAN_SOLAR_AZIMUTH_ANGLE'))
                    .set('solar_zenith',image.get('MEAN_SOLAR_ZENITH_ANGLE') )
                    .set('system:time_start',image.get('system:time_start'));
}

// Function to get clear-sky pixels from Sentinel-2 image
function s2ClearSky(image) {
  var scl = image.select('SCL').multiply(10000);
  var clearSkyPixels = scl.eq(4).or(scl.eq(5)).or(scl.eq(6)).or(scl.eq(11));
  return image.updateMask(clearSkyPixels).divide(10000).copyProperties(image, ["system:time_start"]);
}
// Function to add various vegetation indices to Sentinel-2 imagery
var addVariables = function(image) {
  var DOY = image.date().getRelative('day', 'year');
  var year = image.date().get('year');
  return image
    // Add a NDVI band.
    .addBands(image.normalizedDifference(['nir', 'red']).toDouble().rename('NDVI'))
    // Add a EVI band.
    .addBands(image.expression('2.5*((nir-red)/(nir+6*red-7.5*blue+1))', {
      'nir':image.select('nir'),
      'red':image.select('red'),
      'blue':image.select('blue')
    }).toDouble().rename('EVI'))
    // Add a GCVI
    .addBands(image.expression('nir/green-1',{
      'nir': image.select('nir'),
      'green': image.select('green'),
    }).toDouble().rename('GCVI'))
    // Add a MSAVI2
    .addBands(image.expression('1/2 * (2*nir + 1 - ((2*nir+1)**2 - 8*(nir-red))**(1/2))',{
      'nir': image.select('nir'),
      'red': image.select('red'),
    }).toDouble().rename('MSAVI2'))  
    
    // Add a LSWI band.
    .addBands(image.normalizedDifference(['nir','swir1']).toDouble().rename('LSWI'))
    // Add a NDWI band.
    .addBands(image.normalizedDifference(['green','nir']).toDouble().rename('NDWI'))
    // Add a NDSI band.
    .addBands(image.normalizedDifference(['green','swir1']).toDouble().rename('NDSI'))
    
    // Add NDSVI
    .addBands(image.normalizedDifference(['swir1','red']).toDouble().rename('NDSVI'))
    // Add NDTI
    .addBands(image.normalizedDifference(['swir1','swir2']).toDouble().rename('NDTI'))
    
    // Add S2 red-edge indices (Sen2-Agri)
    // RENDVI = (nir-red2)/(nir+red2)
    // REP = {705+35*(0.5*(red3+red)-red1)/(red2-red1)}/1000
    // PSRI = (red-blue)/red1
    // CRE = red1/nir
    .addBands(image.normalizedDifference(['nir','red2']).toDouble().rename('RENDVI'))
    
    .addBands(image.expression('(705+35*(0.5*(red3+red)-red1)/(red2-red1))/1000',{
      'red3' : image.select('red3'),
      'red2' : image.select('red2'),
      'red1' : image.select('red1'),
      'red' : image.select('red'),
    }).toDouble().rename('REP'))
    
    .addBands(image.expression('(red-blue)/red1',{
      'red': image.select('red'),
      'red1': image.select('red1'), 
      'blue': image.select('blue'), 
    }).toDouble().rename('PSRI'))
    
    .addBands(image.expression('red1/nir',{
      'red1': image.select('red1'),
      'nir': image.select('nir'),
    }).toDouble().rename('CRE'))

    // add a doy band.
    .addBands(ee.Image(DOY).rename('DOY').toDouble())
    // add a year band.
    .addBands(ee.Image(year).rename('Year').toDouble())
    
    .set('DOY',DOY);
};
// Function to add time bands to Sentinel-2 imagery
function addTimeBands(image) {
  var tstamp = ee.Date(image.get('system:time_start'));
  var tdelta = tstamp.difference(startDay, 'year');
  var img_fitting = image
    .addBands(1)
    .addBands(tdelta.multiply(3*Math.PI).sin())
    .addBands(tdelta.multiply(3*Math.PI).cos())
    .addBands(tdelta.multiply(6*Math.PI).sin())
    .addBands(tdelta.multiply(6*Math.PI).cos())
    .addBands(image.select('NDVI','EVI','LSWI'))
    .toDouble();
  return img_fitting;
}
// Function to merge bands
function mergeBands(image,pre){
  return ee.Image(pre).addBands(image);
}


// Crop map generation
// Take Baekdu Mountain Zone (BM) in 2022 as an example
var s2Tools = require("users/nanshany1993/common:sentinel2");
var bands = ee.List(['red2','swir1','swir2','NDVI','EVI','LSWI','NDSVI','NDTI','RENDVI','REP','constant','constant_1','constant_2','constant_3','constant_4']); // 10 index
var aoi = geometry;
var region = geometry;
var year = 2021;

// Construct feature candicates from Sentinel-2 SR images
var startDay = ee.Date.fromYMD(year,1,1);
var endDay = ee.Date.fromYMD(year+1,1,1);
var s2 = ee.ImageCollection("COPERNICUS/S2_SR")
   .filterBounds(aoi)
   .filterDate(startDay,endDay)
   .map(sentinel2sr)
   .map(s2ClearSky)
   .map(addVariables) 
   .map(addTimeBands);
   
// 10-day composite
var startDoy = startDay.getRelative('day','year');
var endDoy = endDay.advance(-1,'day').getRelative('day','year');
var starts = ee.List.sequence(startDoy, endDoy, 10);
var composites = ee.ImageCollection(starts.map(function(start) {
  var doy = start;
  var filtered = s2filtered.filter(ee.Filter.dayOfYear(start, ee.Number(start).add(10))).median().clip(region);
  var bandLength = filtered.bandNames().length();
  var mask = ee.Algorithms.If({                   // mask must be done for time band
    condition : ee.Number(bandLength).gt(0),
    trueCase : filtered.select(0).mask(),
    falseCase : ee.Image(0).clip(region)    
  });
  return filtered.addBands(ee.Image.constant(doy).rename('doy').float())
                .updateMask(mask)
                .set('system:time_start',ee.Date.fromYMD(year,1,1).advance(doy,'day').millis())
                .set('doy',doy)
                .set('length',bandLength)   ;
  }));
  
// Linear interpolation
var size = composites.size();
var LIC = composites.toList(size);
var interpolated = ee.ImageCollection(ee.List.sequence(12,27,1).map(function(i){
  var i = ee.Number(i)
  var before = ee.ImageCollection.fromImages(LIC.slice(i.subtract(9),i))
    .filter(ee.Filter.gt('length',0)).mosaic();
  var after = ee.ImageCollection.fromImages(LIC.slice(i.add(1),i.add(10)).reverse())
    .filter(ee.Filter.gt('length',0)).mosaic();
  var boforeY = before.select(bands);
  var beforedoy = before.select('doy');
  var afterY = after.select(bands);
  var afterdoy = after.select('doy');
  var targetImg = ee.Image(LIC.get(i));
  var currentdoy = ee.Image.constant(targetImg.get('doy')).float();
  var Y = afterY.subtract(boforeY).divide(afterdoy.subtract(beforedoy))
      .multiply(currentdoy.subtract(beforedoy)).add(boforeY);
  var filledImage = ee.Image(ee.Algorithms.If({
    condition : ee.Number(targetImg.get('length')).gt(0), 
    trueCase : targetImg.select(bands).unmask(Y),
    falseCase : Y
  }));
  return filledImage.unmask(0).clip(region)
    .set('system:time_start',targetImg.get('system:time_start'),'doy',targetImg.get('doy')); // can not simply copy all properties of composites
}))  ;

// SG smoothing
var window_size = 7;
var order = 3;
var sgs = s2Tools.sgsmooth(interpolated,bands, order, window_size);

// NDVI/LSWI maximum composite
var NDVIs = sgs.qualityMosaic('NDVI');
var LSWIs = sgs.qualityMosaic('LSWI');
var day10s = ee.Image(sgs.iterate(mergeBands, ee.Image([]))).addBands(NDVIs).addBands(LSWIs);

// Harmonic regression
// y = a + b1*cos(3pi*t) + b2*sin(3pi*t) + b3*cons(6pi*t) +b4*sin(6pi*t) 
var dependent = ee.List(['NDVI','EVI','LSWI']);
var harmonicIndependents = ee.List(['constant', 'constant_2', 'constant_1', 'constant_4' , 'constant_3']);
// The output of the regression reduction is a [X,Y] array image.
var harmonic = s2
  .select(harmonicIndependents.cat(dependent))
  .reduce(ee.Reducer.linearRegression(harmonicIndependents.length(), dependent.length()));
var coefficients = harmonic.select('coefficients').matrixTranspose()
  .arrayFlatten([dependent,harmonicIndependents]).clip(region);  
  
// Dry land mask
var finalImage = day10s.addBands(coefficients).clip(region);
var finalImage = finalImage.addBands(image.select('b1').rename('B1'));
var landmask = finalImage.select('B1').eq(12);
var finalImage = finalImage.updateMask(landmask);

// Test points sample regions
var test_points = ee.FeatureCollection(table2);
var test_points = test_points.randomColumn('random'); 
// print(test_points.size());
// print(test_points.limit(10));
var test_points = test_points.map(function(feat){
  return ee.Feature(feat.geometry(), { 
    grid_code: feat.get('classification'),
  });
});
var test_points = finalImage.sampleRegions({
  collection:test_points,
  properties:['grid_code'],
  scale:10,
  geometries:true,
  tileScale:16
});

// Stratified sampling
var cornSamples = test_points.filter(ee.Filter.eq('grid_code',1));
var nonCornSamples = test_points.filter(ee.Filter.eq('grid_code',2));
var seed = 1;
var cornSamplesStratified = cornSamples.randomColumn('stratified', seed);
var nonCornSamplesStratified = nonCornSamples.randomColumn('stratified', seed);
var numSamples = 500;
var numSamplesfortrain = numSamples* 0.8;
var numSamplesforvalidation = numSamples* 0.2;
var cornTraining = cornSamplesStratified.filter(ee.Filter.lt('stratified', 0.8)).randomColumn('random');
// .limit(numSamplesfortrain);
var nonCornTraining = nonCornSamplesStratified.filter(ee.Filter.lt('stratified', 0.8)).randomColumn('random');
// .limit(numSamplesfortrain);
var cornValidation = cornSamplesStratified.filter(ee.Filter.gte('stratified', 0.8)).randomColumn('random');
// .limit(numSamplesforvalidation);
var nonCornValidation = nonCornSamplesStratified.filter(ee.Filter.gte('stratified', 0.8)).randomColumn('random');
// .limit(numSamplesforvalidation);

var trainingSamples = cornTraining.merge(nonCornTraining);
// print(trainingSamples.size())
var validationSamples = cornValidation.merge(nonCornValidation);
// print(validationSamples.size())
// Export.table.toAsset({
//   collection: validationSamples,
//   description: 'vali_points03'
// })
// Export.table.toAsset({
//   collection: trainingSamples,
//   description: 'train_points03'
// })

// TrainTables is the training dataset, which contains the crop type and Sentinel-2 features of each crop training sample
var trainingSamples = trainingSamples.filterBounds(aoi);
var cropFeatures = ['REP_11', 'REP_15', 'REP_17', 'REP_9', 'REP_10', 'REP_12',  'REP_8', 'NDTI_10', 'swir1_12', 'RENDVI_17', 
'RENDVI_10', 'LSWI_9', 'LSWI_constant_3', 'EVI_9', 'swir2_8', 'swir2_17', 'NDSVI_6','LSWI_8','swir1_8',
'LSWI_16','swir1_9'];

// Train the classifier
var classifier = ee.Classifier.smileRandomForest({
  variablesPerSplit:10,
  numberOfTrees: 100,  
  bagFraction: 0.8
}).train({
  features: trainingSamples,
  classProperty: 'grid_code',
  inputProperties: cropFeatures
});
// print('importance',rf.explain());
var featureImportance = classifier.explain();
// print('importance',featureImportance);
var featureCollection = ee.FeatureCollection([
  ee.Feature(null, featureImportance)
]);

// Export feature importance to Google Drive
Export.table.toDrive({
  collection: featureCollection,
  description: 'FeatureImportance_BM',
  folder: 'FeatureImportance',
  fileNamePrefix: 'FeatureImportance_BM',
  fileFormat: 'CSV'
});

// Classify the image with the trained classifier
var classified = finalImage.classify(classifier);
var test = validationSamples.classify(classifier);
var confusionMatrix = test.errorMatrix('grid_code', 'classification');
var exportAccuracy = ee.Feature(null, {matrix: confusionMatrix.array()});

// Export results
Export.table.toDrive({
  collection: ee.FeatureCollection(exportAccuracy),
  description: 'Accuracy_BM',
  folder:'matrix',
  fileFormat: 'CSV'
});
Export.image.toDrive({
  image: classified,
  // description: label,
  description: 'cropclass_BM',
  folder:'2021',
  scale: 10,
  region: region,  
  maxPixels : 1e13
});
