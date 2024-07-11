
# This script reads the temperature log from a csv file and writes a fortran module that contains the temperature data.
# The module contains a subroutine that takes a time as input and returns the temperature at that time with interpolation.


fileName = 'Reactor_TPcheck.csv'
TemperatureColumn = 1
TimeColumn = 0
startline = 0
outputFileName = 'Temperature.f90'

f = open(fileName, 'r+')
fnew = open(outputFileName, 'w+')

lines = f.readlines()

Times = []
Temperatures = []

for i in range(len(lines)):
    if i < startline:
        continue
    else:
        line = lines[i].split(',')
        try:
            Temperatures.append(float(line[TemperatureColumn]))
            Times.append(float(line[TimeColumn]))
        except:
            if len(Temperatures) > len(Times):
                Temperatures.pop(-1)
            elif len(Temperatures) < len(Times):
                Times.pop(-1)

f.close()

timelines = " &["
Temperaturelines = " &["

for i in range(len(Times)):
    timelines += str(Times[i]) + "_WP, "
    if (len(timelines.split('\n')[-1]) > 90):
        timelines += "&\n     &"
timelines = timelines[:-2]
timelines += "]"

for i in range(len(Temperatures)):
    Temperaturelines += str(Temperatures[i]) + "_WP, "
    if (len(Temperaturelines.split('\n')[-1]) > 90):
        Temperaturelines += "&\n     &"
Temperaturelines = Temperaturelines[:-2]
Temperaturelines += "]"

text = """
module temperature_mod
    use precision, only: WP
    implicit none

    integer, parameter :: npoints = {0}

    real(WP), dimension(npoints) :: times = &
    {1}
    real(WP), dimension(npoints) :: Ts = &
    {2}
    
contains

    subroutine get_temperature(T, time)
        implicit none
        real(WP), intent(out) :: T
        real(WP), intent(in) :: time

        integer :: index
        real(WP) :: frac

        ! use binary search to find the time index
        do index = 1, size(times)-1
            if (time >= times(index) .and. time < times(index+1)) exit
        end do

        frac = (time - times(index)) / (times(index+1) - times(index))
        T = Ts(index) + frac * (Ts(index+1) - Ts(index))


    end subroutine get_temperature
    
end module temperature_mod

""".format(len(Times), timelines, Temperaturelines)

fnew.write(text)
fnew.close()