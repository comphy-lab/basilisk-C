# script-version: 2.0
import sys
print("Catalyst Python Backend: Initializing script...", file=sys.stderr)
# CRITICAL: Force ParaView to initialize in Headless OSMesa Software mode!
sys.argv.append('--osmesa')
sys.argv.append('--force-offscreen-rendering')
from paraview.simple import *

# 1. Initialize Catalyst Options
from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'TimeStep'
options.EnableCatalystLive = 1
options.CatalystLiveTrigger = 'TimeStep'

# 2. Extract Data from Memory (matches the "grid" channel in C)
extractgrid = PVTrivialProducer(registrationName="grid")

# 3. Setup Render View and Visualization
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1920, 1080]
renderView1.Background = [0.1, 0.1, 0.15] # Dark gray background
renderView1.AxesGrid = 'GridAxes3DActor'

# Display the mesh, colored by our scalar 'p'
extractgridDisplay = Show(extractgrid, renderView1, 'UnstructuredGridRepresentation')
extractgridDisplay.Representation = 'Surface'
extractgridDisplay.ColorArrayName = ['CELLS', 'p']

# Center camera automatically
ResetCamera(renderView1)
Render(renderView1)

# 4. Optional: Save Images (Uncomment to enable)
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
pNG1.Trigger = 'TimeStep'
pNG1.Writer.FileName = 'RenderView_{timestep:06d}.png'
pNG1.Writer.ImageResolution = [1920, 1080]
pNG1.Writer.Format = 'PNG'

# 5. Execute
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    SaveExtractsUsingCatalystOptions(options)
