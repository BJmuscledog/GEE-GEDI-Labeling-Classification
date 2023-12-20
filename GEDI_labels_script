 /********************************************************************************************
Take Baekdu Mountain Zone (BM) in 2021 as an example
*********************************************************************************************/ 

// Define image and table variables
var image1 = ee.Image("landuse_2020"); // 30m*30m land cover data (http://www.resdc.cn/)
var image2 = ee.Image("CDL_NE_2017_DYY"); // 30m*30m You et al.(2021) crop type dataset in Northeast China(https://www.geodata.cn/)
var image3 = ee.Image("CDL_NE_2018_DYY"); // 30m*30m You et al.(2021) crop type dataset in Northeast China(https://www.geodata.cn/)
var image4 = ee.Image("CDL_NE_2019_DYY"); // 30m*30m You et al.(2021) crop type dataset in Northeast China(https://www.geodata.cn/)
var image5 = ee.Image("CDL_NE_2020_DYY"); // 30m*30m You et al.(2021) crop type dataset in Northeast China(https://www.geodata.cn/)
var image6 = ee.Image("CDL_NE_2021_DYY"); // 30m*30m You et al.(2021) crop type dataset in Northeast China(https://www.geodata.cn/)
var table1 = ee.FeatureCollection("China_dongbei"); // Administrative boundaries of agricultural areas in Northeast China
var table2 = ee.FeatureCollection("decision-tree-BM"); // Trained decision trees

// Set the region of interest
var BM = ee.Geometry.Point([127.5477075829258, 42.91374280135512]);
var geometry = table1.filterBounds(BM);

// Clip image and apply masks
var img = image6.addBands(image1.select('b1').rename('B1')).updateMask(image1.select('b1').eq(12));
var img = img.addBands(image5.select('b1').rename('b2'))
          .addBands(image4.select('b1').rename('b3'));
var landmask1 = img.select('b2').eq(image5);
var landmask2 = img.select('b3').eq(image4);
var img = img.updateMask(landmask1).updateMask(landmask2);

// GEDI L2A quality mask function
var conditionMask_L2A = function(image){
  var img = img.select('rh0',	'rh1',	'rh2',	'rh3',	'rh4',	'rh5',	'rh6',	'rh7',	'rh8',	'rh9',	'rh10',	'rh11',	'rh12',	'rh13',	'rh14',	'rh15',	'rh16',	'rh17',	
  'rh18',	'rh19',	'rh20',	'rh21',	'rh22',	'rh23',	'rh24',	'rh25',	'rh26',	'rh27',	'rh28',	'rh29',	'rh30',	'rh31',	'rh32',	'rh33',	'rh34',	'rh35',	'rh36',	'rh37',	'rh38',	
  'rh39',	'rh40',	'rh41',	'rh42',	'rh43',	'rh44',	'rh45',	'rh46',	'rh47',	'rh48',	'rh49',	'rh50',	'rh51',	'rh52',	'rh53',	'rh54',	'rh55',	'rh56',	'rh57',	'rh58',	'rh59',
  'rh60',	'rh61',	'rh62',	'rh63',	'rh64',	'rh65',	'rh66',	'rh67',	'rh68',	'rh69',	'rh70',	'rh71',	'rh72',	'rh73',	'rh74',	'rh75',	'rh76',	'rh77',	'rh78',	'rh79',	'rh80',	
  'rh81',	'rh82',	'rh83',	'rh84',	'rh85',	'rh86',	'rh87',	'rh88',	'rh89',	'rh90',	'rh91',	'rh92',	'rh93',	'rh94',	'rh95',	'rh96',	'rh97',	'rh98',	'rh99',	'rh100',
'num_detectedmodes','elev_lowestmode','digital_elevation_model','elev_highestreturn','digital_elevation_model_srtm','beam');
  var mask1 = img.updateMask(img.select('rh100').lt(10).and(img.select('rh100').gte(0)));
  var mask2 = mask1.updateMask(img.select('num_detectedmodes').neq(0));
  var mask3 = mask2.updateMask((img.select('elev_lowestmode').subtract(img.select('digital_elevation_model')).abs().lte(100)));
  var mask4 = mask3.updateMask((img.select('elev_lowestmode').subtract(img.select('digital_elevation_model_srtm')).abs().lte(100)));
  var mask5 = mask4.updateMask(img.select('beam').eq(5).or(img.select('beam').eq(6)).or(img.select('beam').eq(8)).or(img.select('beam').eq(11)));
  return mask5;
};
var qualityMask_L2A = function(img) {
  return img.updateMask(img.select('quality_flag').eq(1))
      .updateMask(img.select('degrade_flag').eq(0));
};

// Create GEDI L2A image collection
var dataset_L2A = ee.ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY')
  .map(qualityMask_L2A)
  .filterDate('2021-7-01', '2021-9-30')
  .filterBounds(geometry)
  .map(conditionMask_L2A)
  .map(function(img) {
    return img.addBands([
      img.expression("(rh100 - rh70)", {'rh100': img.select('rh100'), 'rh70': img.select('rh70')}).rename(['vhd']),
      img.expression("(rh100 - rh0)", {'rh100': img.select('rh100'), 'rh0': img.select('rh0')}).rename(['ghd']),
      img.expression("(rh20 - rh0)", {'rh20': img.select('rh20'), 'rh0': img.select('rh0')}).rename(['rh020'])
    ]);
  })
  .map(function(img) {
    return img.updateMask(img.select('B1').eq(12));
  });
// Display GEDI L2A on the map
// var gediVis = {
//   min: 1,
//   max: 60,
//   palette: 'darkred,red,orange,green,darkgreen',
// };
// Map.addLayer(dataset_L2A.select('beam'), gediVis, 'gedi_l2a');

// GEDI L2B quality mask function
var qualityMask_L2B = function(image) {
  return image.select('l2b_quality_flag').eq(1)
    .and(image.select('degrade_flag').eq(0))
    .and(image.select('algorithmrun_flag').eq(1));
};

// Create GEDI L2B image collection
var dataset_L2B = ee.ImageCollection('LARSE/GEDI/GEDI02_B_002_MONTHLY')
  .map(qualityMask_L2B)
  .filterDate('2021-7-01', '2021-9-30')
  .filterBounds(geometry)
  .map(landmask_l2b)
  .map(function(img) {
    return img.updateMask(img.select('B1').eq(12));
  });
var landMask_L2b = function(image){
  var img0 = image.addBands(img2.select('B1','b1'));
  var img1 = img0.addBands(img3.select('num_detectedmodes','elev_lowestmode','digital_elevation_model','elev_highestreturn','digital_elevation_model_srtm','rh100'));
  var mask0 = img1.updateMask(img1.select('B1').eq(12));
  var mask1 = mask0.updateMask(img1.select('rh100').lt(10).and(img1.select('rh100').gte(0)));
  var mask2 = mask1.updateMask(img1.select('num_detectedmodes').neq(0));
  var mask3 = mask2.updateMask((img1.select('elev_lowestmode').subtract(img1.select('digital_elevation_model')).abs().lte(100)));
  var mask4 = mask3.updateMask((img1.select('elev_lowestmode').subtract(img1.select('digital_elevation_model_srtm')).abs().lte(100)));
  return mask4;
};
var dataset_L2B = dataset_L2B.map(landMask_L2b); 
var img_L2B = dataset_L2B.mosaic().clip(geometry);

// Merge the GEDI L2A and L2B bands
var dataset = dataset_L2A.map(function(img) {
  return img.addBands(img_L2B.select('fhd_normal'));
});

// Set the null value to a proxy value to filter the crushing footprints
var dataset_filter = dataset.map(function(img) {
  return img.unmask(5);
});
var img_filter = dataset_filter.mosaic();
var img = dataset.mosaic();

// // Display GEDI L2B on the map
// Map.addLayer(dataset.select('rh98'), gediVis, 'gedi_l2b');

// Create labels
var img_band = ee.Image.constant(0).rename('new') ;
var img = img.addBands(img_band.select('new'));
var landmask_types_others = img.select('b1').neq(1);
var landmask_types_corn = img.select('b1').eq(1);
var labels_others = img.updateMask(landmask_types_others).clip(geometry).stratifiedSample({
  region: geometry,
  classBand: 'new',
  scale: 25,
  numPoints: 1000, // specified number of labels
  seed: 1,
  dropNulls: true,
  tileScale: 16,
  geometry: true
});

var labels_corn = img.updateMask(landmask_types_corn).clip(geometry).stratifiedSample({
  region: geometry,
  classBand: 'new',
  scale: 25,
  numPoints: 1000, // specified number of labels
  seed: 1,
  dropNulls: true,
  tileScale: 16,
  geometry: true
});
var train_labels = labels_others.merge(labels_corn);

// Select bands for training
var bands = ['rh0', 'rh10', 'rh20', 'rh30', 'rh70', 'rh80', 'rh90', 'rh100', 'vhd', 'ghd', 'rh5', 'rh95', 'fhd_normal'];
var img = img.select(bands);

// Create training labels
var ft = img.sampleRegions({
  collection: train_labels,
  properties: ['b1'],
  scale: 25,
  geometry: true,
  tileScale: 16
});

// Set buffer condition and filter
var pro_labels = ft.map(function(point) {
  var buffer = point.buffer(58).geometry();
  var Clipped = img_filter.clip(buffer);
  var Mean = Clipped.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: buffer,
    scale: 30
  });
  var result = Mean.get('b1');
  var newPoint = point.set('mean', result);
  return newPoint;
}).map(function(feat) {
  return ee.Feature(feat).set('b1', ee.Algorithms.If(ee.Number(feat.get('b1')).eq(1), 1, 0));
}).filter(ee.Filter.eq('mean', 1));

// Randomly shuffle points
var withRandom = pro_labels.randomColumn('random');
// Random Forest Classifier
var rf = ee.Classifier.smileRandomForest({
  numberOfTrees: 100,
  bagFraction: 0.8
}).train({
  features: withRandom,
  classProperty: 'b1',
  inputProperties: bands
});

// // Print feature importance
// print('bands_importance', rf.explain());

// // Apply Random Forest classification
// var point_classification = withRandom.classify(rf);

// // Encode and export decision trees
// var trees = ee.List(ee.Dictionary(rf.explain()).get('trees'));
// var dummy = ee.Feature(geometry, {});
// var col = ee.FeatureCollection(trees.map(function(x) {
//   return dummy.set('tree', x);
// }));
// Export.table.toAsset({
//   collection: col,
//   description: 'decision-trees-BM',
//   assetId: 'users/your_username/decision-trees-BM'
// });

// Decode decision trees
var trained =  ee.FeatureCollection(table2);
var decisionTrees = decode(trained);
var classifier = ee.Classifier.decisionTreeEnsemble(decisionTrees);
function decode(featureCollection) {
  return featureCollection.map(function(feature) {
    var dict = feature.toDictionary();
    var keys = dict.keys().map(function(key) {
      return ee.Number.parse(ee.String(key));
    });
    var value = dict.values().sort(keys).join();
    return ee.Feature(null, { value: value });
  }).aggregate_array('value').join().decodeJSON();
}

// Apply the trained decision tree to the current year label
var classifier = rf.setOutputMode('RAW');
var result = img.classify(classifier);
var prob = result.arrayReduce(ee.Reducer.mean(), ee.List([0])).arrayGet(0);
// print('classifier', prob);
var point_classification = point_classification.map(function(feat) {
  return ee.Feature(feat.geometry(), { grid_code: feat.get('classification') });
});

// Filter labels with low confidence
var conditional = function(feat) {
  return ee.Algorithms.If(
    ee.Number(feat.get('classification')).gte(1).and(ee.Number(feat.get('classification')).lte(1.2)).or(
      ee.Number(feat.get('classification')).gte(1.8)).and(ee.Number(feat.get('classification')).lte(2)),
    feat.set({ classification: 1 }),
    feat
  );
};
var point_classification = point_classification.map(conditional);
var point_classification = point_classification.filter(ee.Filter.eq('classification', 1));

// Export classified points to Google Drive
Export.table.toDrive({
  collection: point_classification,
  description: 'BM_labels_2021',
  folder: 'Labels_2021'
});
