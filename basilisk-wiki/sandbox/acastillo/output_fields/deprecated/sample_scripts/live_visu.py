# script-version: 2.0
# Catalyst state generated using paraview version 5.10.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1920, 1080]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.Background = [1.0, 1.0, 1.0]
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.0, 0.0, 0.0]
renderView1.KeyLightIntensity = 0.65
renderView1.KeyLightElevation = 0.0
renderView1.KeyLightAzimuth = 0.0
renderView1.MaintainLuminance = 1
renderView1.UseAmbientOcclusion = 1
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [1.2982168370389915, 0.09160149132667272, 3.0825954378610856]
renderView1.CameraViewUp = [0.021085850857911898, 0.9990335132975866, -0.03856716481515019]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.8
#renderView1.EnableRayTracing = 1
#renderView1.BackEnd = 'OSPRay raycaster'
renderView1.Shadows = 1
#renderView1.AmbientSamples = 20
#renderView1.SamplesPerPixel = 20
#renderView1.ProgressivePasses = 5
#renderView1.LightScale = 0.8
#renderView1.TemporalCacheSize = 0


# get the material library
materialLibrary1 = GetMaterialLibrary()
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitle = ''
renderView1.AxesGrid.YTitle = ''
renderView1.AxesGrid.ZTitle = ''
renderView1.AxesGrid.XAxisUseCustomLabels = 1
renderView1.AxesGrid.YAxisUseCustomLabels = 1
renderView1.AxesGrid.ZAxisUseCustomLabels = 1

SetActiveView(None)


# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1920, 1080)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)

# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVTrivialProducer'
extractgrid = PVTrivialProducer(registrationName="grid")

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------
# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------


# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from extractgrid
extractgridDisplay = Show(extractgrid, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'pressure'
pressureLUT = GetColorTransferFunction('pressure')
pressureLUT.RGBPoints = [1.448999493606095, 0.231373, 0.298039, 0.752941, 683.439789073575, 0.865003, 0.865003, 0.865003, 1365.430578653544, 0.705882, 0.0156863, 0.14902]
pressureLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'pressure'
pressurePWF = GetOpacityTransferFunction('pressure')
pressurePWF.Points = [0.0, 0.0, 0.5, 0.0, 100, 1.0, 0.5, 0.0]
pressurePWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
extractgridDisplay.Representation = 'Volume'
extractgridDisplay.ColorArrayName = ['CELLS', 'pressure']
extractgridDisplay.LookupTable = pressureLUT
extractgridDisplay.SelectTCoordArray = 'None'
extractgridDisplay.SelectNormalArray = 'None'
extractgridDisplay.SelectTangentArray = 'None'
extractgridDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
extractgridDisplay.SelectOrientationVectors = 'None'
extractgridDisplay.ScaleFactor = 0.1
extractgridDisplay.SelectScaleArray = 'None'
extractgridDisplay.GlyphType = 'Arrow'
extractgridDisplay.GlyphTableIndexArray = 'None'
extractgridDisplay.GaussianRadius = 0.005
extractgridDisplay.SetScaleArray = [None, '']
extractgridDisplay.ScaleTransferFunction = 'PiecewiseFunction'
extractgridDisplay.OpacityArray = [None, '']
extractgridDisplay.OpacityTransferFunction = 'PiecewiseFunction'
extractgridDisplay.DataAxesGrid = 'GridAxesRepresentation'
extractgridDisplay.PolarAxes = 'PolarAxesRepresentation'
extractgridDisplay.ScalarOpacityFunction = pressurePWF
extractgridDisplay.ScalarOpacityUnitDistance = 0.04626252703571701
extractgridDisplay.OpacityArrayName = ['CELLS', 'pressure']
extractgridDisplay.SelectMapper = 'Resample To Image'

# setup the color legend parameters for each legend in this view

# get color legend/bar for pressureLUT in view renderView1
pressureLUTColorBar = GetScalarBar(pressureLUT, renderView1)
pressureLUTColorBar.Title = 'pressure'
pressureLUTColorBar.ComponentTitle = ''

# set color bar visibility
pressureLUTColorBar.Visibility = 1

# show color legend
extractgridDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
pNG1.Trigger = 'TimeStep'

# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'RenderView1_{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = [1920, 1080]
pNG1.Writer.Format = 'PNG'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(pNG1)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'TimeStep'
options.EnableCatalystLive = 1
options.CatalystLiveTrigger = 'TimeStep'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
