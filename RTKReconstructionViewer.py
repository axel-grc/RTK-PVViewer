# Paraview import
from paraview.simple import *
from paraview import vtk
# RTK import
from itk import RTK as rtk
# Common import
import sys
import numpy as np
from math import atan2

# Helper function returning unit vector
def unit_vector(vector):
    """
    Returns the unit vector of the input vector.
    """
    return vector / np.linalg.norm(vector)

def GetCornerPoints(dataBounds, p, u, v):
    """
    Compute the 4 corners point position.

    :param dataBounds: Input data bounds
    :param p: Input data position
    :param u: Input data u direction
    :param v: Input data v direction

    :return list of 3D points
    """

    # Compute half extent in the U and V direction
    eU = 0.5 *(dataBounds[1] - dataBounds[0])
    eV = 0.5 *(dataBounds[3] - dataBounds[2])

    # Compute corner points
    p11 = [p[0] + u[0] * eU + v[0] * eV, p[1] + u[1] * eU + v[1] * eV, p[2] + u[2] * eU + v[2] * eV] # [1,1]
    p10 = [p[0] + u[0] * eU - v[0] * eV, p[1] + u[1] * eU - v[1] * eV, p[2] + u[2] * eU - v[2] * eV] # [1,-1]
    p01 = [p[0] - u[0] * eU + v[0] * eV, p[1] - u[1] * eU + v[1] * eV, p[2] - u[2] * eU + v[2] * eV] # [-1,1]
    p00 = [p[0] - u[0] * eU - v[0] * eV, p[1] - u[1] * eU - v[1] * eV, p[2] - u[2] * eU - v[2] * eV] # [-1,-1]

    return [p11, p10, p01, p00]

def UpdateFrame(timeStep):
    """
    Update animation frame

    param: timeStep: Current step in the animation sequence. [0.0-1.0]
    """
    global geometry, sourceDisplay, projectionsDisplay, projections, line0, line1, line2, line3
    
    # Compute current slice
    numberOfProjections = len(geometry.GetGantryAngles())
    currentSlice = round((timeStep * (numberOfProjections - 1)))

    # Update source position
    sourcePosition = geometry.GetSourcePosition(currentSlice)
    sourceDisplay.Position = [sourcePosition[0], sourcePosition[1], sourcePosition[2]]
    
    # Update current slice
    projectionsDisplay.Slice = currentSlice

    # Compute slice position
    m = geometry.GetProjectionCoordinatesToFixedSystemMatrix(currentSlice)
    bounds = projections.GetDataInformation().DataInformation.GetBounds()
    projorigin = np.array([0.5 *(bounds[0]+bounds[1]), 0.5 *(bounds[2]+bounds[3]), 0.5 *(bounds[4]+bounds[5])] )
    projorigin = np.append(projorigin, [1])
    p = m*projorigin
    
    # Compute slice orientation
    projorigin[0] += 1
    u = m*projorigin - p
    projorigin[0] -= 1
    projorigin[1] += 1
    v = m*projorigin - p
    projorigin[1] -= 1
    # Compute Y angle
    dot = vtk.vtkMath.Dot([1.0,0.0,0.0],unit_vector([u[0], u[1], u[2]]))
    cross = [0,0,0]
    vtk.vtkMath.Cross([1.0,0.0,0.0],unit_vector([u[0], u[1], u[2]]),cross)
    det = vtk.vtkMath.Dot([0,1,0],cross)
    rY = angle = atan2(det, dot) * 180 / 3.1415
    # Compute X angle
    dot = vtk.vtkMath.Dot([0.0,-1.0,0.0],unit_vector([v[0], v[1], v[2]]))
    cross = [0,0,0]
    vtk.vtkMath.Cross([0.0,-1.0,0.0],unit_vector([v[0], v[1], v[2]]),cross)
    det = vtk.vtkMath.Dot([0,0,1],cross)
    rX = angle = atan2(det, dot) * 180 / 3.1415

    # Update slice orientation
    projectionsDisplay.Orientation = [rX,-rY,0]

    # Source-Detector frustum
    projectionCorners = GetCornerPoints(bounds, p, u, v)

    line0.Point1 = sourceDisplay.Position
    line0.Point2 = projectionCorners[0]

    line1.Point1 = sourceDisplay.Position
    line1.Point2 = projectionCorners[1]

    line2.Point1 = sourceDisplay.Position
    line2.Point2 = projectionCorners[2]

    line3.Point1 = sourceDisplay.Position
    line3.Point2 = projectionCorners[3]

    # Update slice position
    offset = -currentSlice * (bounds[5]-bounds[4]) / (numberOfProjections-1) - bounds[4] # WARNING: + bounds[5] vs - bounds[4]
    sourceDetectorDirection = [p[0] - sourcePosition[0], p[1] - sourcePosition[1], p[2] - sourcePosition[2]]
    sourceDetectorDirection = unit_vector(sourceDetectorDirection)
    projectionsDisplay.Position = [p[0] + offset * sourceDetectorDirection[0], p[1] + offset * sourceDetectorDirection[1], p[2] + offset * sourceDetectorDirection[2]]

def callback(caller, *args):
    """
    Callback called on AnimationCueTickEvent.
    Compute the current time step and update animation frame.
    WARNING: The computed time step should be in between 0 and 1

    param: caller: Paraview's callback caller parameter
	param: args: Paraview's callback arguments
    """
    global animationScene

    timeStep = (animationScene.AnimationTime - animationScene.StartTime) / (animationScene.EndTime - animationScene.StartTime)

    if GetTimeTrack().Enabled:
        timeStep = (animationScene.TimeKeeper.Time - animationScene.StartTime) / (animationScene.EndTime - animationScene.StartTime)

    UpdateFrame(timeStep)

from tkinter import Tk
from tkinter.filedialog import askopenfilename
def OpenFile(fileTypes = "", title = ""):
    """
    Open file dialog and return selected filename.

    param: fileTypes: Type of file to look for
	param: title: Dialog title

	return: Selected filename string
    """
    root = Tk()
    root.withdraw()
    fileName = askopenfilename(
        filetypes = fileTypes,
        title = title)
    root.update()
    return fileName

# Main
if not hasattr(sys, 'argv'):
    sys.argv  = ['']

if len ( sys.argv ) > 1:
    print(sys.argv[1])

# Open input files using tkinter
projectionsFileName = OpenFile(
    fileTypes =(("MetaImage", "*.mha *.mhd"),("All Files","*.*")),
    title = "Open Projections file ...")

geometryFileName = OpenFile(
    fileTypes =(("Geometry File", "*.xml"),("All Files","*.*")),
    title = "Open Geometry file ...")

volumeFileName = OpenFile(
    fileTypes =(("MetaImage", "*.mha *.mhd"),("All Files","*.*")),
    title = "Open Volume file ...")

# Disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Render view
renderView = GetActiveViewOrCreate('RenderView')
renderView.ViewSize = [970, 485]

## Projections ##
projections = MetaFileSeriesReader(FileNames=[projectionsFileName], guiName='Projections')
projectionsDisplay = Show(projections, renderView)
projectionsDisplay.SetRepresentationType('Slice')

## Geometry ##
geometryReader = rtk.ThreeDCircularProjectionGeometryXMLFileReader.New()
geometryReader.SetFilename(geometryFileName)
geometryReader.GenerateOutputInformation()
geometry = geometryReader.GetOutputObject()

## Source ##
source = Sphere(guiName='Source')
sourceDisplay = Show(source, renderView)
sourceDisplay.Scale =[50,50,50]

## Volume ##
volumeUY = MetaFileSeriesReader(FileNames=[volumeFileName], guiName='Volume_XY')
volumeDisplayXY = Show(volumeUY, renderView)
volumeDisplayXY.SetRepresentationType('Slice') # 'Volume'
volumeDisplayXY.SliceMode = 'XY Plane' # 'YZ Plane' / 'XZ Plane'

volumeVZ = PassArrays(volumeUY,  guiName='Volume_YZ')
volumeDisplayYZ = Show(volumeVZ, renderView)
volumeDisplayYZ.SetRepresentationType('Slice')
volumeDisplayYZ.SliceMode = 'YZ Plane'

volumeUZ = PassArrays(volumeUY,  guiName='Volume_XZ')
volumeDisplayXZ = Show(volumeUZ, renderView)
volumeDisplayXZ.SetRepresentationType('Slice')
volumeDisplayXZ.SliceMode = 'XZ Plane'

## Source-Detector frustum ##
line0 = Line(guiName='Line0')
line1 = Line(guiName='Line1')
line2 = Line(guiName='Line2')
line3 = Line(guiName='Line3')

lineDisplay0 = Show(line0, renderView)
lineDisplay1 = Show(line1, renderView)
lineDisplay2 = Show(line2, renderView)
lineDisplay3 = Show(line3, renderView)

SetActiveSource(projections)

## Coloring ##
colorTransferFunctionRGBPoints = [-0.060333251953125, 0.231373, 0.298039, 0.752941, 1.366178035736084, 0.865003, 0.865003, 0.865003, 2.792689323425293, 0.705882, 0.0156863, 0.14902]
opacityTransferFunctionPoints = [0.0, 1.0, 0.5, 0.0, 0.4430769085884094, 0.6000000238418579, 0.5, 0.0, 0.8369230628013611, 0.0, 0.5, 0.0, 0.9292307496070862, 1.0, 0.5, 0.0, 1.070769190788269, 0.0, 0.5, 0.0, 2.0, 0.987500011920929, 0.5, 0.0]

ColorBy(GetDisplayProperties(volumeUY), 'MetaImage', True)
GetColorTransferFunction('MetaImage',volumeDisplayXY, True).RGBPoints = colorTransferFunctionRGBPoints
GetOpacityTransferFunction('MetaImage',volumeDisplayXY, True).Points = opacityTransferFunctionPoints
volumeDisplayXY.RescaleTransferFunctionToDataRange(False, True)# Rescale color and/or opacity maps used to eUactly match data range
volumeDisplayXY.LookupTable.EnableOpacityMapping = 1
GetColorTransferFunction('MetaImage',volumeDisplayXY, True).ApplyPreset('Grayscale', True)

volumeDisplayYZ.LookupTable = volumeDisplayXY.LookupTable
volumeDisplayXZ.LookupTable = volumeDisplayXY.LookupTable

volumeDisplayYZ.OpacityTransferFunction = volumeDisplayXY.OpacityTransferFunction
volumeDisplayXZ.OpacityTransferFunction = volumeDisplayXY.OpacityTransferFunction

ColorBy(GetDisplayProperties(projections), 'MetaImage', True)
GetColorTransferFunction('MetaImage',projectionsDisplay, True).RGBPoints = colorTransferFunctionRGBPoints
GetOpacityTransferFunction('MetaImage',projectionsDisplay, True).Points = opacityTransferFunctionPoints
projectionsDisplay.RescaleTransferFunctionToDataRange(False, True)# Rescale color and/or opacity maps used to eUactly match data range
projectionsDisplay.LookupTable.EnableOpacityMapping = 1
GetColorTransferFunction('MetaImage',projectionsDisplay, True).ApplyPreset('X Ray', True)

## Animation ##
animationScene = GetAnimationScene()
animationScene.EndTime = 1
animationScene.FramesPerTimestep = 1
animationScene.Loop = 0
animationScene.PlayMode = 'Sequence' # 'Real Time' / 'Sequence'
animationScene.Duration = 7
animationScene.NumberOfFrames = len(geometry.GetGantryAngles())

# Setup animation callback
scene = animationScene.GetClientSideObject()
callbackid = scene.AddObserver("AnimationCueTickEvent", callback)

# Start rendering
Show()
RenderAllViews()

# Init object position
UpdateFrame(0.0)

# Update the view to ensure updated data information
renderView.ResetCamera()
renderView.Update()

# Clip bounds hack
rv = GetRenderView()
rv.MaxClipBounds = [-5000, 5000, -5000, 5000, -5000, 5000]
cam = GetActiveCamera()
cr = cam.GetClippingRange()
cam.SetClippingRange(cr[0], cr[0] + 5000)
rv.LockBounds = 1
Render()

# Start animation
animationScene.Play()

# Allow for interacting with the scene when using pvpython
Interact()
