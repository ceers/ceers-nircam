import os
import sys

from snowball_run_pipeline import detector1_with_snowball_correction

#################################
### SET I/O AND OTHER SETUP HERE 
UNCALDIR = 'uncals'
OUTPUTDIR = 'calibrated'
# Fraction of available cores to use for multi-processing for Ramp Fitting 
# step (argument = 'maximum_cores'). Options are 'none' (no multi-processing),
# ‘quarter’, ‘half’, and ‘all’. 
MAXCORES = 'half'
# Set DELINTERIM to True to delete interim files to save space - *trapsfilled, 
# *fitopt (if saved), *rateints (unnecessary, CEERS exposures have 1 int)
DELINTERIM = True
#################################


def run_detector1_and_snowballs(dataset, inputdir, outputdir, maxcores):
    """A wrapper to call Harry's snowball_run_pipeline"""
    
    # Run the script on a single uncal
    detector1_with_snowball_correction(dataset, input_dir=inputdir,
                                       output_dir=outputdir, maxcores=maxcores)

    # Rate file is output as *_0_rampfitstep.fits
    # Rateints file is output as *_1_rampfitstep.fits
    rampfit = os.path.join(outputdir, '%s_0_rampfitstep.fits'%dataset)
    rate = os.path.join(outputdir, '%s_rate.fits'%dataset)

    print('Moving %s to %s'%(rampfit, rate))
    os.rename(rampfit, rate)

    if DELINTERIM:
        # Now delete some of the interim files to save space
        #   trapsfilled - The number of filled traps at each pixel at the
        #                 end of the exposure, output by the Persistence step
        #   fitopt - Optional output with fit info for each ramp segment 
        #            and pixel
        #   rampfitstep - The rateints file that contains slope, etc., for each
        #                 integration. CEERS exposures have only 1 integration
        trapsfilled = os.path.join(outputdir, '%s_trapsfilled.fits'%dataset)
        fitopt = os.path.join(outputdir, '%s_fitopt.fits'%dataset)
        rampfitstep = os.path.join(outputdir, '%s_1_rampfitstep.fits'%dataset)
        for f in [trapsfilled, fitopt, rampfitstep]:
            try:
                os.remove(f)
            except OSError:
                print('File %s doesnt exist'%f)


def main():
    if len(sys.argv) != 2:
        print('Provide name of uncal file')
        exit()
    uncalfile = sys.argv[1]

    # Get dataset ID string for interim file outputs
    dataset = uncalfile.split('_uncal.fits')[0]
    run_detector1_and_snowballs(dataset, UNCALDIR, OUTPUTDIR, MAXCORES)
                


if __name__ == '__main__':
    main()

