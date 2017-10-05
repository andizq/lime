import os

try:
  from taskinit import *
except ImportError:
  class DummyCasaLog:
    originStr=''

    def origin(self, originStr):
      self.originStr = originStr

    def post(self, errStr, priority='Message'):
      print "%s: %s" % (priority, errStr)

  casalog = DummyCasaLog()

def raytrace(gridInFile,moldatfile,dust,filename,imgres,pxls,unit,freq,rotationStyle
  ,theta,phi,incl,posang,azimuth,distance,nThreads,traceRayAlgorithm,doLine,nchan
  ,velres,trans,molI,bandwidth,source_vel,doInterpolateVels,polarization):

  # Note: rotationStyle and doLine are defined at this level purely to define some logical constraints on the subset of parameters listed; tey are not passed directly to LIME.

  casalog.origin('raytrace')

  import threading
  import time

  try:
    import modellib as ml
  except ImportError as err:
    casalog.post("Cannot find module 'modellib'", priority='ERROR')
    return

  try:
    import lime
  except ImportError as err:
    casalog.post("Cannot find module 'lime'", priority='ERROR')
    return

  class LimeExistentialLog:
    complete   = False
    errorState = False
    errorMessage = ''

  class LimeThread (threading.Thread):
    def __init__(self, par, images, limeLog):
      self._par = par
      self._images = images
      self._limeLog = limeLog
      threading.Thread.__init__(self)

    def run(self):
      try:
        lime.runLime(self._par, self._images)
      except Exception, e:
        self._limeLog.errorState = True
        self._limeLog.errorMessage = e

      self._limeLog.complete = True

  sleepSec = 2#10
  gridDoneMessageIsPrinted = False
  mainDoneMessageIsPrinted = False
  raysDoneMessageIsPrinted = False
  nItersCounted = 0

  fnDict = {'dust':dust,'gridInFile':gridInFile}
  for variableName in fnDict.keys():
    fileName = fnDict[variableName]
    if not fileName is None and fileName!='':
      if not os.path.exists(fileName):
        casalog.post("Cannot find %s file %s" % (variableName, fileName), priority='ERROR')
        return

  if not isinstance(moldatfile, list):
    moldatfile = [moldatfile]

  for i in range(len(moldatfile)):
    fileName = moldatfile[i]
    if not fileName is None and fileName!='':
      if not os.path.exists(fileName):
        casalog.post("Cannot find moldatfile[%d] file %s" % (i, fileName), priority='ERROR')
        return

  #vvvvvvvvvvvvvvvvvvvv
  # Bodging up a dummy modellib setting, since Lime won't run without it.
  #
  ml.setCurrentModel('CG97')

  # Set all parameters:
  ml.setParamDouble('Mstar', 0.5)     # Msun
  ml.setParamDouble('Rstar', 2.0)     # Rsun
  ml.setParamDouble('Tstar', 4000.)   # K
  ml.setParamDouble('bgdens', 1.e-4)  # 1/cm^3
  ml.setParamDouble('hph', 4.0)
  ml.setParamDouble('plsig1', -1.0)
  ml.setParamDouble('rin', 1.0)       # AU
  ml.setParamDouble('rout', 100.0)    # AU
  ml.setParamDouble('sig0', 0.01)     # g/cm^2

  # Model CG97 does not provide the results 'abundance',
  # 'doppler', 'temperature, and 'bmag'. Hence these
  # Results have to be provided by Functions

  # Constant abundance (1.e-9):
  ml.setFunction('abundance', 'scalarConst')
  ml.setFunctionParamDouble('abundance', 'val', 1.e-9)

  # Constant doppler (100 m/s):
  ml.setFunction('doppler', 'scalarConst')
  ml.setFunctionParamDouble('doppler', 'val', 100.)

  # Zero bmag:
  ml.setFunction('bmag', 'vectorConstR')
  ml.setFunctionParamDouble('bmag', 'val', 0.)

  # Set temperature identical to t_dust:
  ml.setTempIdentTdust()

  # Done setting up Model Library:
  ml.finalizeConfiguration()
  if not ml.isFinalizedConfiguration():
    casalog.post("Error in finalizing the model config.", priority='ERROR')
    return

  #^^^^^^^^^^^^^^^^^^^^^

  # Set input parameters for lime:
  par = lime.createInputPars()
#  par.radius            = 2000.0*macros["AU"]
#  par.minScale          = 0.5*macros["AU"]
#  par.pIntensity        = 4000
#  par.sinkPoints        = 3000
  par.dust              = dust
#  par.outputfile        = "populations.pop"
#  par.binoutputfile     = "restart.pop"
#  par.gridfile          = "grid.vtk"
#  par.pregrid           = "pregrid.asc"
#  par.restart           = "restart.pop"
  par.gridInFile        = gridInFile
  par.collPartIds        = [1] # must be a list, even when there is only 1 item.
  par.nMolWeights        = [1.0] # must be a list, even when there is only 1 item.
#  par.collPartNames     = ["phlogiston"] # must be a list, even when there is only 1 item.
#  par.collPartMolWeights = [2.0159] # must be a list, even when there is only 1 item.
#  par.gridDensMaxValues = [1.0] # must be a list, even when there is only 1 item.
#  par.gridDensMaxLoc    = [[0.0,0.0,0.0]] # must be a list, each element of which is also a list with 3 entries (1 for each spatial coordinate).
#  par.tcmb              = 2.72548
#  par.lte_only          = False
#  par.init_lte          = False
#  par.samplingAlgorithm = 0
#  par.sampling          = 2 # Now only accessed if par.samplingAlgorithm==0 (the default).
#  par.blend             = False
  par.polarization      = polarization
  par.nThreads          = nThreads
  par.nSolveIters       = 0
  par.traceRayAlgorithm = traceRayAlgorithm
#  par.resetRNG          = False
  par.doSolveRTE        = False
#  par.gridOutFiles      = ['','','','',"grid_5.ds"]
  par.moldatfile        = moldatfile[:]
#  par.girdatfile        = ["myGIRs.dat"] # must be a list, even when there is only 1 item.

  # Define the images:
  images = []

  # Image #0:
  img = lime.createImage()
  img.nchan             = nchan
  img.trans             = trans
  img.molI              = molI
  img.velres            = velres
  img.imgres            = imgres
  img.pxls              = pxls
  img.unit              = unit
  img.freq              = freq
  img.bandwidth         = bandwidth
  img.source_vel        = source_vel
  img.theta             = theta
  img.phi               = phi
  img.incl              = incl
  img.posang            = posang
  img.azimuth           = azimuth
  img.distance          = distance
  img.doInterpolateVels = doInterpolateVels
  img.filename          = filename
#  img.units             = "0,1"

  images.append(img)

  limeLog = LimeExistentialLog()

  lime.setSilent(False)

  # Run LIME in its own thread:
  limeThread = LimeThread(par, images, limeLog)
  limeThread.start()

  casalog.post('Starting LIME run.')

  # Loop with delay - every so often, check LIME status and print it to the casa logger.
  while not limeLog.complete:
    time.sleep(sleepSec)
    limeStatus = lime.getStatus()

    if limeStatus.error!=0:
      break

    if limeStatus.statusGlobal!=0:
      casalog.post('LIME run is complete')
      break

    if limeStatus.statusRayTracing==0:
      if nItersCounted<=0:
        if limeStatus.statusGridBuilding==0:
          pass
        elif not gridDoneMessageIsPrinted:
          gridDoneMessageIsPrinted = True
          casalog.post('Grid is complete')

      if nItersCounted<par.nSolveIters:
        while nItersCounted<=limeStatus.numberIterations:
          nItersCounted += 1
          casalog.post('Iteration %d/%d' % (nItersCounted, par.nSolveIters))
          casalog.post('Min SNR %e  median %e' % (limeStatus.minsnr, limeStatus.median))

    elif not raysDoneMessageIsPrinted:
      raysDoneMessageIsPrinted = True
      casalog.post('Raytracing is complete')

  limeStatus = lime.getStatus()
  if limeStatus.error!=0:
    casalog.post(limeStatus.message, priority='ERROR')

  elif limeLog.errorState:
    casalog.post(limeLog.errorMessage, priority='ERROR')

  # If we get to here, LIME is finished or has crashed, so return.
  return



