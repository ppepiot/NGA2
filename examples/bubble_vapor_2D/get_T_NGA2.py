# Import libraries
import os
from paraview.simple import *

desired_times = [0.5, 0.6, 0.7, 0.8]
tolerance = 0.01
res_path = '/Users/shayanhbi/Repositories/NGA2/examples/bubble_vapor_2D/temperature'
if not os.path.exists(res_path):
    try:
        os.mkdir(res_path)
    except OSError as error:
        print(f"Error: {error}")

# Find source
source = FindSource('nga.case')

# Get the TimeKeeper object
timeKeeper = GetTimeKeeper()

# Get the number of time steps
num_time_steps = len(timeKeeper.TimestepValues)

# Get the field
T = source.CellData.GetArray('Temperature')

# Index of the plot over line
line_ind = -1

# Iterate over each time step
for time_step in range(num_time_steps):

    # Set the current time step
    SetActiveView(GetRenderView())
    animationScene = GetAnimationScene()
    animationScene.AnimationTime = timeKeeper.TimestepValues[time_step]

    # Get the current time
    time = animationScene.AnimationTime
    print('')

    # Check if the current time is one of the specified times
    if any(abs(time - desired_time) <= tolerance for desired_time in desired_times):

        # Increment the index
        line_ind = line_ind + 1

        # create a new 'Plot Over Line'
        center_line = PlotOverLine(registrationName=f"{'ctr_'}{line_ind}", Input=source)
        center_line.SamplingPattern = 'Sample At Segment Centers'
        center_line.Point1 = [0.0, 0.0, 0.0]
        center_line.Point2 = [0.565, 0.0, 0.0]

        # Get the VTK object associated with the plot over line
        vtk_object = center_line.GetClientSideObject()

        # Force the update of the pipeline
        vtk_object.Update()

        # Get the VTK output data object from the VTK object
        output_data = vtk_object.GetOutput()

        # Get the CellData
        point_data = output_data.GetPointData()
        
        # Get the arrays from CellData
        T_ctr = point_data.GetArray('Temperature')
        x_ctr = point_data.GetArray('arc_length')

        # Write to file
        num_points = T_ctr.GetNumberOfTuples()
        with open(res_path + '/num' + str(desired_times[line_ind]) + '.dat', 'w') as file_ctr:
            file_ctr.write(f"x            T\n")
            for i in range(1, num_points):
                T_val = T_ctr.GetValue(i)
                if T_val > 0.0: 
                    file_ctr.write(f"{x_ctr.GetValue(i)}            {T_ctr.GetValue(i)}\n")
        file_ctr.close()

        # Print message
        print('Extracted data for t = '+str(time))

    # Delet the line
    Delete(center_line)
    del center_line