###
###
###
# THIS NEED TO BE CHANGED
# this has the folder that contains the merged images per well
# can put as many as you would like - it will go over one by one
# it will also skip any images it has already processed
folder_paths = [
  "/vast/scratch/users/moore.z/projects/incu/data/20231107/data/form/composite/", 
  "/vast/scratch/users/moore.z/projects/incu/data/20231107/data/test/composite/"
  ]
###
###
###

import sys

from os.path import expanduser

import glob, os

import java.io.File as File

from ij import IJ
from ij import WindowManager

from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.io import TmXmlWriter
from fiji.plugin.trackmate.io import CSVExporter
from fiji.plugin.trackmate.detection import LogDetectorFactory
from fiji.plugin.trackmate.cellpose import CellposeDetectorFactory
from fiji.plugin.trackmate.cellpose import CellposeSettings
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
# from fiji.plugin.trackmate.tracking import LAPUtils
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter

# we have to do the following to avoid errors with utf8 char
reload(sys)
sys.setdefaultencoding('utf-8')

def run_pipeline(image_file):
  
  well = os.path.splitext(os.path.basename(image_file))[0]
  
  # this needs to be full path
  imp = IJ.openImage(image_file)
  imp.show()
  
  #----------------------------
  # create the model
  #----------------------------
  
  model = Model()
  
  # send all messages to log
  model.setLogger(Logger.IJ_LOGGER)
  
  #------------------------
  # prepare settings object
  #------------------------
  
  settings = Settings(imp)
  
  # configure detector - We use the Strings for the keys
  settings.detectorFactory = CellposeDetectorFactory()
  settings.detectorSettings['TARGET_CHANNEL'] = 1
  settings.detectorSettings['OPTIONAL_CHANNEL_2'] = 0
  settings.detectorSettings['CELLPOSE_MODEL'] = CellposeSettings.PretrainedModel.CUSTOM
  
  ###
  ###
  ###
  # THESE NEED TO BE CHANGED
  # this tells the script where your environment is as well as the model
  settings.detectorSettings['CELLPOSE_PYTHON_FILEPATH'] = '/vast/scratch/users/moore.z/envs/cellpose_segment/bin/python'
  settings.detectorSettings['CELLPOSE_MODEL_FILEPATH'] = '/vast/scratch/users/moore.z/projects/incu/process/2_segmentation/model/CP_20230323_110216'
  ###
  ###
  ###
  
  settings.detectorSettings['CELL_DIAMETER'] = 100.0
  settings.detectorSettings['USE_GPU'] = True
  settings.detectorSettings['SIMPLIFY_CONTOURS'] = True
  
  # configure spot filters
  filter1 = FeatureFilter('QUALITY', 30, True)
  settings.addSpotFilter(filter1)
  
  # Configure tracker
  settings.trackerFactory = SparseLAPTrackerFactory()
  # settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap()
  settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
  settings.trackerSettings['LINKING_MAX_DISTANCE' ] 		= 200.0
  settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE' ]	= 200.0
  settings.trackerSettings['MAX_FRAME_GAP']	= 2
  settings.trackerSettings['ALTERNATIVE_LINKING_COST_FACTOR'] = 1.05
  settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = False
  settings.trackerSettings['ALLOW_TRACK_MERGING'] = True
  settings.trackerSettings['SPLITTING_MAX_DISTANCE'] = 15.0
  
  # add ALL the feature analyzers
  settings.addAllAnalyzers()
  
  #-------------------
  # instantiate plugin
  #-------------------
  
  trackmate = TrackMate(model, settings)
  
  #--------
  # process
  #--------
  
  ok = trackmate.checkInput()
  if not ok:
      print("not ok input")
      return
  #    sys.exit(str(trackmate.getErrorMessage()))
  
  ok = trackmate.process()
  if not ok:
      print("not ok process")
      return
  #    sys.exit(str(trackmate.getErrorMessage()))
  
  #----------------
  # results
  #----------------
  
  # selection.
  selectionModel = SelectionModel( model )
  
  # read the default display settings.
  ds = DisplaySettingsIO.readUserDefault()
   
  # echo results with the logger we set at start:
  # model.getLogger().log( str( model ) )
  
  #----------------
  # export results
  #----------------
  
  # save xml
  out_file_xml = File(folder + "/" + well + ".xml")
  writer = TmXmlWriter(out_file_xml) 
  writer.appendSettings(settings)		      # settings content
  writer.appendModel(model)			          # tracking results
  writer.appendDisplaySettings(ds)	      # display settings
  writer.writeToFile()
  
  # export spot csv
  out_file_csv = folder + "/" + well + "-spots.csv"
  CSVExporter.exportSpots(out_file_csv, model, True)
  print("done!")

for folder in folder_paths:
    print("path is " + folder)
    
    tiff_paths = [
        os.path.join(folder, f) 
        for f in os.listdir(folder) 
        if f.endswith(('.tif', '.tiff'))
    ]
    
    for tiff in tiff_paths:
        spot_file = os.path.splitext(tiff)[0] + "-spots.csv"
        
        if os.path.exists(spot_file):
            print("calculation done for " + tiff)
            continue
        
        if tiff.endswith((".tif", ".tiff")):
            print("running " + tiff)
            run_pipeline(image_file = tiff)
            
            os.system("find /tmp/TrackMate* -mmin +5 -delete")
            
            
sys.exit()
